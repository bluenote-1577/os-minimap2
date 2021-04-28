#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/*
 * Inserts a k-mer representing the window around a syncmer 
 * @param t     1-indexed position offset t for syncmer
 * 
 *
 *
 */
static inline void insert_syncmer (int len, const char *str, int t, uint32_t rid, int i, int s, int k, mm128_v *p, void *km, int strand){
    int start_pos = i - k + 1;
    int end_pos = start_pos + k;
    if (end_pos > len){
        fprintf(stderr,"end pos, len, start_pos,i = %d,%d,%d,%d",end_pos,len,start_pos,i);
        //assert(end_pos < len);
        return;
    }
    if (start_pos < 0){
        fprintf(stderr,"end pos, len, start_pos, i = %d,%d,%d,%d",end_pos,len,start_pos,i);
        assert(start_pos > -1);
    }
    mm128_t info = { UINT64_MAX, UINT64_MAX };
    int is_ambiguous = 0;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
    for (int j = start_pos; j < end_pos; ++j){
		int c = seq_nt4_table[(uint8_t)str[j]];
        if (c >= 4){
            is_ambiguous = 1;
            break;
        }
        kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
    }
    if (kmer[0] == kmer[1]) return; // skip "symmetric k-mers" as we don't know it strand
    if (is_ambiguous) return;
    //int z = kmer[0] < kmer[1]? 0 : 1; // strand

    info.x = hash64(kmer[strand], mask) << 8 | k;
    info.y = (uint64_t)rid<<32 | (uint32_t)(end_pos-1)<<1 | strand;
    kv_push(mm128_t, km, *p, info);

}

/**
 * Find symmetric (k,s,t)-open syncmers on a DNA sequence. Note: only takes syncmer
 * if the last minimum in the window is at position t. can modify such that
 * takes syncmer if any minimium is at position t. 
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param s      s-mer size
 * @param t      t offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void open_sync_sketch(void *km, const char *str, int len, int k, int s, int t, uint32_t rid, int is_hpc, mm128_v *p)
{
    uint64_t shift1 = 2 * (s - 1), mask = (1ULL<<2*s) - 1, smer[2] = {0,0};
    uint64_t kmer[2] = {0,0}, kmer_mask = (1ULL<<2*k) - 1, shift_kmer = 2*(k-1);
	int i, j, l, buf_pos, min_pos, smer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;
    int w = k - s + 1;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
    assert(k > s);
    //Don't want to deal with hpc right now
    assert(!is_hpc);
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

    //Take k-mer from pos i
    //look into how buffer_pos interacts with t. Itl ooks like a circular buffer sort of thing for a linear window
    //try to understand what /smer_span l is doing. i think l keeps track of the k-mer that's building, and when it's
    //>= s then we hash the k-mer. After it's greater than s, I think it just becomes some variable hwich keeps increasing
    //which becomes 0 afer a bad character, signaling the start of the process again. 
	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
        int strand = 0;
		if (c < 4) { // not an ambiguous base
			int z;
            //Not sure if HPC works for syncmers
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				smer_span += skip_len;
				if (tq.count > s) smer_span -= tq_shift(&tq);
			} else smer_span = l + 1 < s? l + 1 : s;
			smer[0] = (smer[0] << 2 | c) & mask;           // forward s-mer
			kmer[0] = (kmer[0] << 2 | c) & kmer_mask;           // forward k-mer
			smer[1] = (smer[1] >> 2) | (3ULL^c) << shift1; // reverse s-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift_kmer; // reverse k-mer
			if (smer[0] == smer[1]) continue; // skip "symmetric s-mers" as we don't know it strand
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = smer[0] < smer[1]? 0 : 1; // strand
            strand = kmer[0] < kmer[1]? 0 : 1;
			++l;
			if (l >= s && smer_span < 256) {
				info.x = hash64(smer[z], mask) << 8 | smer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, smer_span = 0;
        //buf_pos is the last s-mer added. buf_pos+1 is the first k-mer in the window. 
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        //fprintf(stderr,"info %d,%d, s = %d, t = %d, smer = %s \n",info.x,info.y,s,t);
        int offset = (strand == 0)? t : w - t + 1;
        //int offset = t;
		if (info.x <= min.x) { // a new minimum; then write the old min
			//if (l >= w + s && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
            if (offset == w && l >= w + s - 1 && min.x != UINT64_MAX){
                insert_syncmer(len,str,offset,rid,i,s,k,p,km,strand);
            };

		} 
        else if (buf_pos == min_pos) { // old min has moved outside the window
            //fprintf(stderr,"buf_pos,t,w,min_pos = %d,%d,%d,%d",buf_pos,t,w,min_pos);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical s-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest s-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + s - 1 && min.x != UINT64_MAX) { // write identical s-mers
                if ((buf_pos + offset) % w == min_pos){
                    //fprintf(stderr,"inserting syncmer");
                    insert_syncmer(len,str,offset,rid,i,s,k,p,km,strand);
                }
			}
		}
        else if ((buf_pos + offset) % w == min_pos && l >= w+s-1){
            //fprintf(stderr,"inserting syncmer");
            insert_syncmer(len,str,offset,rid,i,s,k,p,km,min.y%2);
        }
		if (++buf_pos == w) buf_pos = 0;
	}
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}
