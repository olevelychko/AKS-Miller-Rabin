# AKS primality test and Miller-Rabin algorithm
This Java project contains the AKS test (pseudocode taken from “PRIMES is in P” By Manindra Agrawal, Neeraj Kayal, and Nitin Saxena) and
Miller-Rabin algorithm (pseudocode taken from “Primality Testing for Beginners” by Lasse Rempe-Gillen and Rebecca Waldecker). This project was created
as a part of my Master’s Thesis, so it is more an example than an efficient realisation. The idea was to compare a deterministic test (AKS) with a random one (Miller-Rabin).

## Requirements
* Java Development Kit (JDK) 17 or later.
* IntelliJ IDEA or any other IDE for Java.

## Running the Program
1. Clone or download this repository.
2. Open the project folder in IntelliJ IDEA.
3. Make sure the project SDK is set to your installed JDK version.
4. Open the Java file located in the src directory.
5. Run the “main()” method.

## Implementation details
SecureRandom() in the main function is used to generate a random number with a specified bit length. You can choose a desired length x by “int bitLength = x;”.

### Miller-Rabin
Miller-Rabin algorithm can be launched by stating “MillerRabin(randomNumb);”; you can substitute “randomNumb” with your number “n”. Note that both Miller-Rabin and AKS
use BigInteger type, so to give your number “n”, you need to write “BigInteger n = new BigInteger(String val, int radix);”. You can also specify the desired number of 
independent rounds in Miller-Rabin by changing the value “k” inside the Miller-Rabin() function (k=1 by default; error probability of the algorithm is (1/4)^k).

### AKS
AKS algorithm can be launched by stating “AKS(randomNumb);”. It can use two different functions for polynomial multiplication: multiplyModPoly() (just usual multiplication
with modular reduction, quite slow) and binarySegmentation() (more efficient method, but can work only with number size less than 156 bits). Switching between these two functions
happens in “selectingMultiplyMethod()”, but you can choose one specific function if you wish.
