import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.*;


public class PrimalityTesting {

    private static BigInteger randomBigInt(BigInteger max) {
        SecureRandom random = new SecureRandom();
        BigInteger number;
        int nBits = max.bitLength();
        do {
            number = new BigInteger(nBits, random);
        } while (number.compareTo(max) >= 0 || number.signum() <= 0);
        return number;
    }


    private static boolean isPerfectPower(BigInteger n) {
        int t = n.bitLength();
        for (int k = 2; k <= t; k++) {
            BigInteger lowBound = BigInteger.TWO;
            BigInteger uppBound = n;
            while (lowBound.compareTo(uppBound) <= 0) {
                BigInteger half = lowBound.add(uppBound).divide(BigInteger.TWO);
                BigInteger powerK = half.pow(k);
                int comparison = powerK.compareTo(n);
                if (comparison == 0 && half.pow(k).equals(n)) return true;
                if (comparison < 0) lowBound = half.add(BigInteger.ONE);
                else uppBound = half.subtract(BigInteger.ONE);
            }
        }
        return false;
    }


    private static BigInteger findR(BigInteger n) {
        double log2n = n.bitLength() - 1;
        double twoPowLog = Math.pow(log2n, 2);
        double threePowLog = Math.pow(log2n, 3);
        BigInteger maxR = BigInteger.valueOf((long) Math.ceil(threePowLog));
        BigInteger maxLog = BigInteger.valueOf((long) Math.ceil(twoPowLog));
        if (maxR.compareTo(BigInteger.valueOf(3)) < 0) return BigInteger.valueOf(3);
        boolean goodOrd;
        for (BigInteger r = BigInteger.TWO; r.compareTo(maxR) < 1; r = r.add(BigInteger.ONE)) {
            if ((r.gcd(n).equals(BigInteger.ONE))) {
                goodOrd = true;
                for (BigInteger i = BigInteger.ONE; i.compareTo(maxLog) < 1; i = i.add(BigInteger.ONE)) {
                    if (n.modPow(i, r).equals(BigInteger.ONE)) {
                        goodOrd = false;
                        break;
                    }
                }
                if (goodOrd) return r;
            }
        }
        return BigInteger.ZERO;
    }


    private static HashMap<BigInteger, BigInteger> multiplyModPoly(HashMap<BigInteger, BigInteger> polyA, HashMap<BigInteger, BigInteger> polyB, BigInteger r, BigInteger n) {
        HashMap<BigInteger, BigInteger> prod = new HashMap<>();
        for (var f : polyA.entrySet()) {
            for (var g : polyB.entrySet()) {
                BigInteger deg = f.getKey().add(g.getKey()).mod(r);
                BigInteger coeff = f.getValue().multiply(g.getValue()).mod(n);
                if (!coeff.equals(BigInteger.ZERO)) {
                    prod.put(deg, prod.getOrDefault(deg, BigInteger.ZERO).add(coeff).mod(n));
                }
            }
        }
        return prod;
    }


    private static HashMap<BigInteger, BigInteger> selectingMultiplyMethod(HashMap<BigInteger, BigInteger> polyA, HashMap<BigInteger, BigInteger> polyB, BigInteger r, BigInteger n){
        if(n.bitLength() > 156) return multiplyModPoly(polyA, polyB, r, n);
          else
            return binarySegmentation(polyA, polyB, r, n);
    }


    private static HashMap<BigInteger, BigInteger> polyFastExp(HashMap<BigInteger, BigInteger> base, BigInteger r, BigInteger n) {
        BigInteger exp = n;
        HashMap<BigInteger, BigInteger> result = new HashMap<>();
        result.put(BigInteger.ZERO, BigInteger.ONE);
        HashMap<BigInteger, BigInteger> basePoly = new HashMap<>(base);
        while (exp.signum() > 0) {
            if (exp.testBit(0)) result = selectingMultiplyMethod(result, basePoly, r, n);
            exp = exp.shiftRight(1);
            if (exp.signum() > 0) basePoly = selectingMultiplyMethod(basePoly, basePoly, r, n);
        }
        result.entrySet().removeIf(entry -> entry.getValue().equals(BigInteger.ZERO));
        return result;
    }


    private static HashMap<BigInteger, BigInteger> binarySegmentation(HashMap<BigInteger, BigInteger> polyA, HashMap<BigInteger, BigInteger> polyB, BigInteger r, BigInteger n) {
        BigInteger maxA = polyA.keySet().stream().reduce(BigInteger.ZERO, BigInteger::max).add(BigInteger.ONE);
        BigInteger maxB = polyB.keySet().stream().reduce(BigInteger.ZERO, BigInteger::max).add(BigInteger.ONE);
        BigInteger maxCoefA = polyA.values().stream().map(BigInteger::abs).max(BigInteger::compareTo).orElse(BigInteger.ZERO);
        BigInteger maxCoefB = polyB.values().stream().map(BigInteger::abs).max(BigInteger::compareTo).orElse(BigInteger.ZERO);
        BigInteger bound = (maxA.max(maxB)).multiply(maxCoefA).multiply(maxCoefB);
        int intBound = bound.bitLength() + 1;
        BigInteger A = numbFromPoly(polyA, intBound);
        BigInteger B = numbFromPoly(polyB, intBound);
        BigInteger C = A.multiply(B);
        BigInteger bitMask = BigInteger.ONE.shiftLeft(intBound).subtract(BigInteger.ONE);
        HashMap<BigInteger, BigInteger> c = new HashMap<>();
        BigInteger maxDeg = maxA.add(maxB).subtract(BigInteger.ONE);
        for (BigInteger i = BigInteger.ZERO; i.compareTo(maxDeg) < 0; i = i.add(BigInteger.ONE)) {
            BigInteger shift = i.multiply(BigInteger.valueOf(intBound));
            if (shift.bitLength() > 31)
                throw new IllegalArgumentException("Exponent is too large for binary segmentation");
            BigInteger cCoeff = C.shiftRight(shift.intValueExact()).and(bitMask).mod(n);
            if (!cCoeff.equals(BigInteger.ZERO)) {
                BigInteger exp = i.mod(r);
                c.merge(exp, cCoeff, (aCoeff, bCoeff) -> aCoeff.add(bCoeff).mod(n));
            }
        }
        return c;
    }


    private static BigInteger numbFromPoly(HashMap<BigInteger, BigInteger> x, int intBound) {
        BigInteger X = BigInteger.ZERO;
        for (Map.Entry<BigInteger, BigInteger> entry : x.entrySet()) {
            BigInteger degX = entry.getKey();
            BigInteger shift = degX.multiply(BigInteger.valueOf(intBound));
            if (shift.bitLength() > 31)
                throw new IllegalArgumentException("Exponent is too large for binary segmentation");
            X = X.add(entry.getValue().shiftLeft(shift.intValueExact()));
        }
        return X;
    }


    private static boolean checkCongruence(BigInteger n, BigInteger r) {
        long logN = n.bitLength() - 1;
        HashMap<BigInteger, BigInteger> baseX = new HashMap<>();
        baseX.put(BigInteger.ONE, BigInteger.ONE);
        HashMap<BigInteger, BigInteger> Xrhs = new HashMap<>();
        HashMap<BigInteger, BigInteger> Xlhs;
        Xrhs.put(n.mod(r), BigInteger.ONE);
        for (BigInteger a = BigInteger.ONE; a.compareTo(r.sqrt().multiply(BigInteger.valueOf(logN))) < 0; a = a.add(BigInteger.ONE)) {
            Xrhs.put(BigInteger.ZERO, a); //X^n+a
            baseX.put(BigInteger.ZERO, a); //X+a
            Xlhs = polyFastExp(baseX, r, n); //(X+a)^n
            if (!Xlhs.equals(Xrhs)) return false;
        }
        return true;
    }


    private static boolean AKS(BigInteger n) {
        System.out.println("Agrawal, Kayal, Saxena algorithm");
        System.out.println("n is " + n);
        if (isPerfectPower(n)) {
            System.out.println(n + " is a perfect power => composite (step 1)");
            return false;
        }
        BigInteger r = findR(n);
        System.out.println("r is " + r + " (step 2)");
        for (BigInteger a = BigInteger.TWO; a.compareTo(r) <= 0; a = a.add(BigInteger.ONE)) {
            BigInteger gcdNum = a.gcd(n);
            if (!(gcdNum.equals(BigInteger.ONE)) && gcdNum.compareTo(n) < 0) {
                System.out.println(n + " is composite by gcd with " + a + " (step 3)");
                return false;
            }
        }
        if (n.compareTo(r) <= 0) {
            System.out.println(n + " is prime (step 4) ");
            return true;
        }
        if (!checkCongruence(n, r)) {
            System.out.println(n + " is composite (step 5)");
            return false;
        }
        System.out.println(n + " is prime (end)");
        return true;
    }


    private static boolean MillerRabin(BigInteger n) {
        System.out.println("Miller-Rabin algorithm");
        if (n.remainder(BigInteger.TWO).equals(BigInteger.ZERO) && n.compareTo(BigInteger.TWO) > 0) {
            System.out.println(n + " is composite");
            return false;
        }
        int s = 0;
        int k = 1;
        BigInteger x;
        BigInteger d = n.subtract(BigInteger.ONE);
        while ((d.remainder(BigInteger.TWO)).equals(BigInteger.ZERO) && !d.equals(BigInteger.ZERO)) {
            d = d.divide(BigInteger.TWO);
            s++;
        }
        for (int j = 0; j < k; j++) {
            x = randomBigInt(n);
            if (!(x.gcd(n).equals(BigInteger.ONE))) {
                System.out.println(n + " is composite");
                return false;
            }
            BigInteger x0 = x.modPow(d, n);
            if (x0.equals(BigInteger.ONE) || x0.equals(n.subtract(BigInteger.ONE))) {
                continue;
            }
            int i = 0;
            for (; i < s; i++) {
                x0 = x0.modPow(BigInteger.TWO, n);
                if (x0.equals(n.subtract(BigInteger.ONE))) {
                    break;
                }
            }
            if (i == s) {
                System.out.println(n + " is composite");
                return false;
            }
        }
        System.out.println(n + " is probably prime");
        return true;
    }


    public static void main(String[] args) throws Exception {
       int bitLength = 36;

        SecureRandom random = new SecureRandom();
        BigInteger randomNumb = new BigInteger(bitLength, random);
        randomNumb = randomNumb.setBit(bitLength - 1);
        System.out.println("The random number = " + randomNumb);
        long start0 = System.nanoTime();
        System.out.println(MillerRabin(randomNumb));
        long end0 = System.nanoTime();
        long duration0 = end0 - start0;
        System.out.printf(Locale.US, "Execution time of Miller-Rabin: %.2f ms%n", duration0 / 1_000_000.0);
        System.out.println();
        long start1 = System.nanoTime();
        System.out.println(AKS(randomNumb));
        long end1 = System.nanoTime();
        long duration1 = end1 - start1;
        System.out.printf(Locale.US, "Execution time of AKS: %.2f ms%n", +duration1 / 1_000_000.0);

    }
}











