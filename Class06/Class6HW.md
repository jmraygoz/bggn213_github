# Class06-HW
Jessica PID A15647602

HW Question 6

``` r
# Can you improve this analysis code?

library(bio3d)
```

    Warning: package 'bio3d' was built under R version 4.3.3

``` r
s1 <- read.pdb("4AKE")  # kinase with drug
```

      Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE")  # kinase no drug
```

      Note: Accessing on-line PDB file
       PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y")  # kinase with drug
```

      Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](Class6HW_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](Class6HW_files/figure-commonmark/unnamed-chunk-1-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](Class6HW_files/figure-commonmark/unnamed-chunk-1-3.png)

HW Question 6 Answer:

``` r
library(bio3d) #This will provide the read.pdb file

#Function name
Analyse_protein <- function(Protein){ 
  
  s1 <- read.pdb(Protein) #This is to access the online pdb file data
  s1.chainA <- trim.pdb(s1, chain="A", elety="CA") #Trim and filter the data to retrieve specific data from chain A and elety CA
  s1.b <- s1.chainA$atom$b #Selecting the data from the previously filtered chain A and only elety CA, then selecting the atom's b values 
  plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")  #plot the selected proteins Bfactor vs residues.
  
  return() #this will return and generate the plot
  
   }
```

``` r
Analyse_protein("1E4Y")
```

      Note: Accessing on-line PDB file

    Warning in get.pdb(file, path = tempdir(), verbose = FALSE):
    /var/folders/n9/p7kz9c9s4dg0dqlk47c9nzx40000gn/T//RtmphbSYU8/1E4Y.pdb exists.
    Skipping download

![](Class6HW_files/figure-commonmark/unnamed-chunk-3-1.png)

    NULL

``` r
#This allows for you to add a Protein PDB file to generate a plot that contains the 
```

``` r
Analyse_protein("4AKE") #It works for different proteins
```

      Note: Accessing on-line PDB file

    Warning in get.pdb(file, path = tempdir(), verbose = FALSE):
    /var/folders/n9/p7kz9c9s4dg0dqlk47c9nzx40000gn/T//RtmphbSYU8/4AKE.pdb exists.
    Skipping download

![](Class6HW_files/figure-commonmark/unnamed-chunk-4-1.png)

    NULL

``` r
Analyse_protein("1AKE") #It worrks again...
```

      Note: Accessing on-line PDB file

    Warning in get.pdb(file, path = tempdir(), verbose = FALSE):
    /var/folders/n9/p7kz9c9s4dg0dqlk47c9nzx40000gn/T//RtmphbSYU8/1AKE.pdb exists.
    Skipping download

       PDB has ALT records, taking A only, rm.alt=TRUE

![](Class6HW_files/figure-commonmark/unnamed-chunk-5-1.png)

    NULL
