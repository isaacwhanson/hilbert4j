
package starlig.ht.hilbert4j;

/**
 * Created by isaac on 7/11/15.
 */
public class Main {

    public static void main(String[] args) {
        byte n = 3;
        byte m = 2;
        Hilbert hil = new Hilbert(n, m);
        for (int i = 0; i < hil.getLength(); i++) {
            // i
            prn(String.format("i=%d ", i));
            // p = H-1(i)
            long p = hil.hilbertInv(i);
            prn("p=" + hil.toBinaryString(p));
            // h = H(p)
            //long h = hil.hilbert(p);
            //prn(String.format(" h=%d ", h));
            // p2 = H-1(h)
            //long p2 = hil.hilbertInv(h);
            //prn("p2=" + hil.toBinaryString(p2) + " ");

            prn("\n");
        }
    }

    public static void prn(Object s) {
        System.out.print(s);
    }
}
