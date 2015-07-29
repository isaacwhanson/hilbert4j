
package starlig.ht.hilbert4j;

/**
 * @author Isaac W Hanson
 */
public class Hilbert {

    /**
     * @param x bit-field
     * @return minimum bits to represent value - 1
     */
    public static byte log2(long x) {
        if (x == 0)
            return 0;
        return (byte) (63 - Long.numberOfLeadingZeros(x));
    }

    /**
     * @param i number of bits
     * @return a bit mask of i lsb
     */
    public static long mask(byte i) {
        return ~(~0L << i);
    }

    /**
     * @param x number
     * @return number of trailing set bits (ones)
     */
    public static byte tsb(long x) {
        return (byte) Long.numberOfTrailingZeros(~x);
    }

    /**
     * @param i index
     * @return ith reflected gray code
     */
    public static long gray(long i) {
        return i ^ i >>> 1;
    }

    /**
     * @param gc reflected gray code
     * @return index in set of gray codes
     */
    public static long grayInv(long gc) {
        byte bits = (byte) (log2(gc) + 1);
        long result = gc;
        for (byte shift = 1; shift < bits; shift++)
            result ^= gc >>> shift;
        return result;
    }
    protected byte n;
    protected byte m;

    /**
     * @param n number of dimensions
     * @param m bits per dimension (order)
     */
    public Hilbert(byte n, byte m) {
        this.n = n;
        this.m = m;
    }

    /**
     * @return number of points, 2^(n*m)
     */
    public long getLength() {
        return 1L << n * m;
    }

    /**
     * @param x long to convert
     * @return binary string
     */
    public String toBinaryString(long x) {
        return String.format("%" + Integer.toString(n * m) + "s",
                             Long.toBinaryString(x)).replace(" ", "0");
    }

    /**
     * @param i index on curve
     * @return intra-sub hypercube dimension
     */
    protected int intra(long i) {
        return i == 0 ? 0 : (byte) Long.remainderUnsigned(tsb(i - 1L + (i & 1L)), n);
    }

    /**
     * @param i index on curve
     * @return entry point on sub-hypercube
     */
    protected long entry(long i) {
        return i == 0L ? 0L : gray(i - 1 >>> 1 << 1);
    }

    /**
     * @param x bit-field
     * @param i rotate bits right
     * @return x with n lsb rotated by i right >>>
     */
    protected long rotateRight(long x, byte i) {
        return ((x >>> i | x << n - i) & mask(n)) | (x & ~mask(n));
    }

    /**
     * @param x bit-field
     * @param i rotate bits left
     * @return x with n lsb rotated by i left <<
     */
    protected long rotateLeft(long x, byte i) {
        return ((x << i | x >>> n - i) & mask(n)) | (x & ~mask(n));
    }

    /**
     * @param e entry point
     * @param d intra-sub hypercube dimension
     * @param b index on curve
     * @return rotated entry on sub-hypercube
     */
    protected long transform(long e, byte d, long b) {
        return rotateRight(b ^ e, (byte) (d + 1));
    }

    /**
     * @param e entry point
     * @param d intra-sub hypercube dimension
     * @param b index on curve
     * @return inverse transform
     */
    protected long transformInv(long e, byte d, long b) {
        return transform(rotateRight(e, (byte) (d + 1)), (byte) (n - d - 2), b);
    }

    /**
     * @param x bit-field
     * @param i set offset in units of n
     * @return bits i*n -> (i*n)+n of x
     */
    protected byte read(long x, byte i) {
        return (byte) (x >>> i * n & mask(n));
    }

    /**
     * @param p interlaced point data
     * @return hilbert index of point p
     */
    public long hilbert(long p) {
        long h = 0L;
        long e = 0L;
        byte d = 0;
        for (byte i = (byte) (m - 1); i > -1; i--) {
            long w = grayInv(transform(e, d, read(p, i)));
            h = h << n | w;
            e ^= rotateLeft(entry(w), (byte) (d + 1));
            d = (byte) ((d + intra(w) + 1) % n);
        }
        return h;
    }

    /**
     * @param h hilbert index
     * @return point on curve at h (interlaced)
     */
    public long hilbertInv(long h) {
        long p = 0L;
        long e = 0L;
        byte d = 0;
        for (byte i = (byte) (m - 1); i > -1; i--) {
            long w = read(h, i);
            p |= transformInv(e, d, gray(w)) << i * n;
            e ^= rotateLeft(entry(w), (byte) (d + 1));
            d = (byte) ((d + intra(w) + 1) % n);
        }
        return p;
    }
}
