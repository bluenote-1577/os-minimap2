# os-minimap2 : minimap2 with open syncmer capabilities 

This is a forked version of minimap2 which has the capability of using open syncmers (Edgar 2021, see https://peerj.com/articles/10805/), instead of minimizers for selecting k-mers to be used for seeding. 

It is known that open syncmers perform better than minimizers for k-mer matching tasks under mutation. We show in Shaw and Yu (2021) that open syncmers perform well against other k-mer selection methods in theory and that it improves sensitivity of mapping, especially when reads are relatively divergent from the reference.

## Quickstart
```
git clone https://github.com/bluenote-1577/os-minimap2
cd os-minimap2 && make
## Using open syncmers with (s,t) = (11,3).
./minimap2 -a --syncs 11 --synct 3 test/MT-human.fa test/MT-orang.fa > test.sam 
```
To use syncmers, specify parameters ``--syncs (s) --synct (t)`` where ``t,s`` correspond to parameters in the definition of open syncmers. We require ``k > s``. If these parameters are not specified, minimap2 will use minimizers instead. 

For detailed instructions on how to use minimap2, see https://github.com/lh3/minimap2. 

Currently homo-polymer compressed runs are not supported. Indexing and mapping as separate steps is also not supported.

## Citations

Shaw, Jim and Yu, Yun William. "Theory of local k-mer selection with applications to
long-read alignment" Bioinformatics (2022).
