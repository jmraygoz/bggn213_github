# Class06: R functions
Jessica PID A15647602

``` r
add <- function(x, y = 1){
  x + y
}
```

Can I use it?

``` r
add(1,1)
```

    [1] 2

``` r
add(x = 1, y = 100)
```

    [1] 101

``` r
add(c(100, 1, 100), y = 1) 
```

    [1] 101   2 101

``` r
add(1, 1)
```

    [1] 2

Q. Make a function; `generate_DNA()` that makes a random nucleotide
sequence of any length

``` r
generate_DNA <- function(){

bases <- c("A", "T", "C", "G")
sequence <- sample(bases, size = 10, replace = TRUE)
}
```

That is my working snipet now i can make a function.

``` r
generate_DNA <- function(length){
  bases <- c("A", "T", "C", "G")
  sequence <- sample(bases, size = length, replace = TRUE)
  return(sequence)
}
```

``` r
generate_DNA(15)
```

     [1] "G" "C" "A" "G" "T" "G" "A" "A" "G" "T" "T" "G" "T" "A" "C"

``` r
aa <- unique(bio3d::aa.table$aa1)[1:20]
```

``` r
generate_protein <- function(length){
  amino_acids<- c(aa)
  protein_sequence <- sample(amino_acids, size = length, replace = TRUE)
  protein_sequence <- paste(protein_sequence, collapse = "") #added to remove space
   return(protein_sequence)
  }
```

``` r
generate_protein(10)
```

    [1] "MYHTGVYQQT"

Q. Generate random protein sequences of length 6 to 12.

``` r
answer <- sapply(6:12, generate_protein)
answer
```

    [1] "TDRHPH"       "KLMSRPW"      "PPNKAAWI"     "WIRQPTNIF"    "LLLKRISKVW"  
    [6] "KQYHKWIQGWI"  "VLDRACGIVFPC"

``` r
cat( paste(">id.", 6:12, "\n",  answer, sep = ""), sep = "\n")
```

    >id.6
    TDRHPH
    >id.7
    KLMSRPW
    >id.8
    PPNKAAWI
    >id.9
    WIRQPTNIF
    >id.10
    LLLKRISKVW
    >id.11
    KQYHKWIQGWI
    >id.12
    VLDRACGIVFPC
