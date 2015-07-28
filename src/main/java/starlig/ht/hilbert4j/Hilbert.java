package starlig.ht.hilbert4j;

/**
 * @author Isaac W Hanson
 */
public class Hilbert {

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
     * @param i rotate bits right
     * @return x with n lsb rotated by i right >>>
     */
    public long rotateRight(long x, int i) {
        return (((x >>> i) | (x << (n - i))) & mask(n)) | (x & ~mask(n));
    }

    /**
     * @param x bit-field
     * @param i rotate bits left
     * @return x with n lsb rotated by i left <<
     */
    public long rotateLeft(long x, int i) {
        return (((x << i) | (x >>> (n - i))) & mask(n)) | (x & ~mask(n));
    }

    /**
     * @param x long to convert
     * @return binary string
     */
    public String toBinaryString(long x) {
        String format = "%" + Integer.toString(n * m) + "s";
        return String.format(format, Long.toBinaryString(x)).replace(" ", "0");
    }

    /**
     * @return number of points, 2^(n*m)
     */
    public long getLength() {
        return 1L << (n * m);
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
        return i == 0 ? 0 : (int) Long.remainderUnsigned(tsb(i - 1L + (i & 1L)), n);
    }

    /**
     * @param i index on curve
     * @return entry point on sub-hypercube
     */
    protected long entry(long i) {
        return i == 0L ? 0L : getGrayCode(((i - 1) >>> 1) << 1);
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
     * @param x bit-field
     * @param i set offset in units of n
     * @param j get offset in units of n
     * @return bits (0 -> n-1) + j * n of x, copied to bits (0 -> n-1) + i * n
     */
    protected long extract(long x, int i, int j) {
        long y = 0;
        for (int k = 0; k < n; k++)
            y = setBit(y, k + (i * n), getBit(x, k + (j * n)));
        return y;
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
            long w = getGrayCodeInv(transform(e, d, extract(p, 0, i)));
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
    public long hilbertInv(long h) {
        long e = 0L;
        int d = 0;
        long p = 0L;
        for (int i = m - 1; i > -1; i--) {
            long w = extract(h, 0, i);
            long l = transformInv(e, d, getGrayCode(w));
            p ^= extract(l, i, 0);
            e ^= rotateLeft(entry(w), d + 1);
            d = (d + intra(w) + 1) % n;
        }
        return p;
    }
}
