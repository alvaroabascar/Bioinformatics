Perl Project
===

This is a Perl project aiming to solve the problem explained below.
The report of this project is stored under /doc.

Getting the Topmost Scoring Sequences from Position Weight Matrices
-----------

<a href="http://en.wikipedia.org/wiki/Position_weight_matrix"><strong>Position Weight Matrices</strong></a> are a simple way to model signals appearing on DNA and protein sequences. They summaryze the frequencies of every letter of the nucleotide or amino acid alphabets at a given position of the signal, for instance a splice site or a transcription factor binding site. This means that we can obtain a score, either from the absolute frequencies or by transforming them into log-likelihoods, for a given sequence to determine if it contains the signal pattern or not. Yet, <strong>they can also be used as sequence generators</strong>. The idea is to write a Perl script to <strong>produce all the possible sequences a PWM can match and report them ranked by score</strong>, for the shake of simplicity, we will work with <em>nucleotides</em>. On some cases, it can be unfeasible to produce all the possible sequences; the worst scenario is a PWM for a random sequence, where the score of any of the four nucleotides is 0.25, as the number of posible sequences will be 4n (being n the length of the matrix). For instance a matrix for a nucleotide signal made of 10 nucleotides, n=10, can generate 410 = 1048576 different sequences. As we do not want to fill the disk with too many sequences, we can <strong>set a cut-off, say here M=10000 sequences. However, we still must produce only the M sequences that have the highest scores.</strong> The simplest way to achieve that is to <em>use a secondary PWM where the positions are reordered by the score of the most frequent nucleotide per position, as well as an array with the original positions, taken as scalar values, recording the order of the new PWM. The best approach can be a recursive function that start generating new sequences iterating over the nucleotides per position sorted from higher frequencies to lower ones. The output records should contain two fields, the generated sequence and the corresponding score, both can be calculated simultaneously.</em> Finally, the script should read a text file containing the PWM in <a href="http://meme.nbcr.net/meme/doc/transfac-format.html">TRANSFAC format</a>, having a section with five fields to store in memory (the ones marked in green, for the position and the relative frequencies for A, C, G and T, in that order). Below you can find two matrices that you have to run separately through your program. Discuss the results obtained from them on the report.
<br /><br /><br />
AC  M00034<br />
XX<br />
ID  V$P53_01<br />
XX<br />
DE  tumor suppressor p53<br />
XX<br />
BF  T00671; p53; Species: Homo sapiens.<br />
BF  T01806; p53; Species: Mus musculus.<br />
XX<br />
P0      A      C      G      T<br />
01      4      0     13      0      G<br />
02      5      0     12      0      G<br />
03     15      0      2      0      A<br />
04      0     17      0      0      C<br />
05     17      0      0      0      A<br />
06      0      0      0     17      T<br />
07      0      0     17      0      G<br />
08      0     13      0      4      C<br />
09      0     17      0      0      C<br />
10      0     17      0      0      C<br />
11      0      0     17      0      G<br />
12      0      0     17      0      G<br />
13      2      0     15      0      G<br />
14      0     17      0      0      C<br />
15     17      0      0      0      A<br />
16      0      0      0     17      T<br />
17      0      0     17      0      G<br />
18      0      2      0     15      T<br />
19      0     13      0      4      C<br />
20      0      7      2      7      Y<br />
XX<br />
//<br />
<br />
<br />
AC  M00097<br />
XX<br />
ID  V$PAX6_01<br />
XX<br />
DE  Pax-6<br />
XX<br />
BF  T00681; Pax-6; Species: Mus musculus.<br />
BF  T01122; Pax-6; Species: Homo sapiens.<br />
XX<br />
P0      A      C      G      T<br />
01     15      7      6     10      N<br />
02     21      9      3     10      N<br />
03     10      9     10     18      N<br />
04      8     14      9     16      N<br />
05      3      2      4     38      T<br />
06      2      0      1     44      T<br />
07      3     29      1     14      C<br />
08     40      5      1      1      A<br />
09      3     39      0      5      C<br />
10      1      0     44      2      G<br />
11      1     36      7      2      C<br />
12     23      2      1     21      W<br />
13      1      4      0     42      T<br />
14      2     13     26      3      G<br />
15     40      1      6      0      A<br />
16     14     11     15      7      N<br />
17      2      4      3     37      T<br />
18      1      0     20     25      K<br />
19     13     17      9      4      N<br />
20     14      8      4      6      N<br />
21      4     12      3      9      N<br />
XX<br />
//<br />
