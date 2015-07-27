package starlig.ht.hilbert4j;

/**
 * @author Isaac W Hanson
 */
public class Hilbert {

    protected int n;
    protected int m;

    /**
     * @param n number of dimensions
     * @param m bits per dimension (order)
     */
    public Hilbert(int n, int m) {
        this.n = n;
        this.m = m;
    }

    /**
     * @param x bit-field
     * @param i position
     * @return bit i in x
     */
    public static int getBit(long x, int i) {
        return (int) ((x >>> i) & 1L);
    }

    /**
     * @param x bit-field
     * @param i position
     * @param b value
     * @return x with b at i
     */
    public static long setBit(long x, int i, int b) {
        long bi = 1L << i;
        return b == 1 ? x | bi : x & ~bi;
    }

    /**
     * @param x bit-field
     * @return minimum bits to represent value - 1
     */
    public static int log2(long x) {
        if (x == 0)
            return 0;
        return 63 - Long.numberOfLeadingZeros(x);
    }

    /**
     * @param i number of bits
     * @return a bit mask of i lsb
     */
    public static long mask(int i) {
        return ~(~0L << i);
    }

    /**
     * @param x number
     * @return number of trailing set bits (ones)
     */
    public static int tsb(long x) {
        return Long.numberOfTrailingZeros(~x);
    }

    /**
     * @param x bit-field
     * @param i rotate bits right
     * @param j lsb to rotate
     * @return x with j lsb rotated by i right >>>
     */
    public static long rotateRight(long x, int i, int j) {
        return (((x >>> i) | (x << (j - i))) & mask(j)) | (x & ~mask(j));
    }

    /**
     * @param x bit-field
     * @param i rotate bits left
     * @param j lsb to rotate
     * @return x with j lsb rotated by i left <<
     */
    public static long rotateLeft(long x, int i, int j) {
        return (((x << i) | (x >>> (j - i))) & mask(j)) | (x & ~mask(j));
    }

    /**
     * @param x long to convert
     * @param i number of places
     * @return binary string
     */
    public static String toBinaryString(long x, int i) {
        String format = "%" + Integer.toString(i) + "s";
        return String.format(format, Long.toBinaryString(x)).replace(" ", "0");
    }

    /**
     * @return number of points, 2^(n*m)
     */
    public long getLength() {
        return 1L << (n * m);
    }

    /**
     * @param i number to convert
     * @return binary string with length n*m
     */
    public String toBinaryString(long i) {
        return toBinaryString(i, n * m);
    }

    /**
     * @param i index
     * @return ith reflected gray code
     */
    public long getGrayCode(long i) {
        return i ^ (i >>> 1);
    }

    /**
     * @param gc reflected gray code
     * @return index in set of gray codes
     */
    public long getGrayCodeInv(long gc) {
        int bits = log2(gc) + 1;
        int shift = 1;
        long result = gc;
        while (shift < bits) {
            result ^= (gc >>> shift);
            shift++;
        }
        return result;
    }

    /**
     * @param i index on curve
     * @return intra-sub hypercube dimension
     */
    protected int intra(long i) {
        return i == 0 ? 0 : (int) Long.remainderUnsigned(inter(i - 1L + (i & 1L)), n);
    }

    /**
     * @param i index on curve
     * @return inter-sub hypercube dimension
     */
    protected int inter(long i) {
        return tsb(i);
    }

    /**
     * @param i index on curve
     * @return entry point on sub-hypercube
     */
    protected long entry(long i) {
        return i == 0L ? 0L : getGrayCode(((i - 1) >>> 2) << 1);
    }

    protected long rotateLeft(long x, int i) {
        return rotateLeft(x, i, n);
    }

    protected long rotateRight(long x, int i) {
        return rotateRight(x, i, n);
    }

    /**
     * @param e entry point
     * @param d intra-sub hypercube dimension
     * @param b index on curve
     * @return rotated entry on sub-hypercube
     */
    protected long transform(long e, int d, long b) {
        return rotateRight(b ^ e, d + 1);
    }

    /**
     * @param e entry point
     * @param d intra-sub hypercube dimension
     * @param b index on curve
     * @return inverse transform
     */
    protected long transformInv(long e, int d, long b) {
        return transform(rotateRight(e, d + 1), n - d - 2, b);
    }

    /**
     * @param p interlaced point data
     * @return hilbert index of point p
     */
    public long hilbert(long p) {
        long h = 0L;
        long e = 0L;
        int d = 0;
        for (int i = m - 1; i > -1; i--) {
            long l = 0;
            for (int j = 0; j < n; j++)
                l = setBit(l, j, getBit(p, j + (i * n)));
            long w = getGrayCodeInv(transform(e, d, l));
            h = (h << n) | w;
            e ^= rotateLeft(entry(w), d + 1);
            d = (d + intra(w) + 1) % n;
        }
        return h;
    }

    /**
     * @param h hilbert index
     * @return point on curve at h (interlaced)
     */
    public long hilbertInverse(long h) {
        long e = 0L;
        int d = 0;
        long p = 0L;
        for (int i = m - 1; i > -1; i--) {
            long w = 0L;
            for (int j = 0; j < n; j++)
                w = setBit(w, j, getBit(h, (i * n) + j));
            long l = transformInv(e, d, getGrayCode(w));
            for (int j = 0; j < n; j++)
                p = setBit(p, j + (i * n), getBit(l, j));
            e ^= rotateLeft(entry(w), d + 1);
            d = (d + intra(w) + 1) % n;
        }
        return p;
    }
}
