# epic

diffuse domain chip seq caller based on SICER.

This is a pre-alpha release. Please do aggressively report issues, quirks, complaints and anything that just feels slightly not right to the issue tracker.

I have been using it with great success, but no others have looked at it so you are likely to encounter some issues.

#### License

MIT

#### Quickstart

```
pip install bioepic
# you only need git clone to get the test data
git clone https://github.com/endrebak/epic.git
epic -i control epic/examples/test.bed epic/examples/control.bed
```

#### Why another one?

MACS2 is great for narrow peaks, but epic performs better on diffuse domains. For medium size domains, such as PolII, our tests indicate that both perform about equally well, but epic uses only a fraction of the time.

#### Why not SICER?

SICER is a wonderful algorithm, but advances in the Python data science libraries has made it possible to implement it much more efficiently. Furthermore, SICER was not made to handle the mountains of data we have now; it simply cannot run on very large datasets due to (sensible) restrictions in the original implementation.

#### Credit

Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao and Weiqun Peng for the original SICER. Please consider citing their paper (*in addition* to our pre-print) if you use epic.

Endre Bakken Stovner for the implementation of epic.

#### Why the name epic?

It stands for electronic pic [sic] caller or epigenome cartographer, whichever you prefer. Or perhaps it isn't just another bogus bioinformatics acronym. Hope you find the name fitting.

#### Bug in the original SICER?

I would appreciate it if anyone can send me the output islands with FDR from running the original SICER on the example data they provide. I seem to remember finding a bug in the original, but cannot be bothered to get SICER running again.
