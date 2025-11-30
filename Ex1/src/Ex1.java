/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
  htt*ps://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};
	/* this mine function for calculate the polynom with every length. There is another one example of the same function< but
	prefer use this.
	*/

    private static double calPolynom(double[] polynom, double x) {
        if (polynom.length < 1) return 0;

        double ans = polynom[0];
        for (int i = polynom.length - 1; i > 0; i--) {
            ans += polynom[i] * Math.pow(x, i);
        }
        return ans;
    }

    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     *
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function computes a polynomial representation from a set of 2D points on the polynom.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Note: this function only works for a set of points containing up to 3 points, else returns null.
     *
     * @param xx
     * @param yy
     * @return an array of doubles representing the coefficients of the polynom.
	 Dicribtion: 
	 This function calculate our polynom by three points. In the beginning I build the polynom from the two points if our function
	 is straight line and after from three points if it porabola. This code uses a formula from the link. there is explaing how to find 
	 vertex's polynom. 
	 in the end I get ppolynom as an ans. 
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {
            if (lx == 2) {
                double x1 = xx[0];
                double y1 = yy[0];
                double x2 = xx[1];
                double y2 = yy[1];

                double m = (y1 - y2) / (x1 - x2);

                double b = y1 - (m * x1);

                ans = new double[]{b, m, 0};
            } else {
                double x1 = xx[0];
                double y1 = yy[0];
                double x2 = xx[1];
                double y2 = yy[1];
                double x3 = xx[2];
                double y3 = yy[2];

                double D = x1 * x1 * (x2 - x3) - x2 * x2 * (x1 - x3) + x3 * x3 * (x1 - x2);
                if (D == 0) return null;

                double Da = y1 * (x2 - x3) - y2 * (x1 - x3) + y3 * (x1 - x2);
                double Db = x1 * x1 * (y2 - y3) - x2 * x2 * (y1 - y3) + x3 * x3 * (y1 - y2);
                double Dc = x1 * x1 * (x2 * y3 - x3 * y2) - x2 * x2 * (x1 * y3 - x3 * y1) + x3 * x3 * (x1 * y2 - x2 * y1);

                double a = Da / D;
                double b = Db / D;
                double c = Dc / D;

                ans = new double[]{c, b, a};
            }
        }
        return ans;
    }

    /**
     * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     *
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true iff p1 represents the same polynomial function as p2.
	 This code uses by my first private function for checking if we have two polunoms that equal mathematically.
	 A diffrence of polynoms has to be fewer than epsilon, ether it will return me false. 
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;

        for (int x = 0; x <= 2; x++) {
            double y1 = Ex1.calPolynom(p1, x);
            double y2 = Ex1.calPolynom(p2, x);

            if (Math.abs(y1 - y2) > EPS) {
                ans = false;
                break;
            }
        }

        return ans;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2.0,3.1,-1.2} will be presented as the following String  "-1.2x^2 +3.1x +2.0"
     *
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function:
	 This function give us polynom as string from an array. 
	 First, we check if we have polynom that equals to 0. If yes, it will return me 0.
	 Secondly, if we have only coeff, so it will return us just a plus without any sign 
	 After, if we have coeff and x< but x's degree is 1, so it will return to ans "x".
	 If we have coeff and x with degree< we return to ans "x^" with i, cuz it will be a degree. 
	 in the end, we just return full ans with full polynom.
     */
    public static String poly(double[] poly) {
        String ans = "";
        if (poly.length == 0) {
            ans = "0";
        } else {
            for (int i = poly.length - 1; i >= 0; i--) {
                if (poly[i] == 0) continue;

                if (!ans.isEmpty()) {
                    ans += poly[i] >= 0 ? " +" : " ";
                }

                if (i == 0) {
                    ans += poly[i];
                } else if (i == 1) {
                    ans += poly[i] + "x";
                } else {
                    ans += poly[i] + "x^" + i;
                }
            }
        }
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
     * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     *
     * @param p1  - first polynomial function
     * @param p2  - second polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     */
    
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans; 

        while ((x2 - x1) > eps) {
            double xm = (x1 + x2) / 2;
            double f1 = Ex1.calPolynom(p1, x1) - Ex1.calPolynom(p2, x1);
            double fm = Ex1.calPolynom(p1, xm) - Ex1.calPolynom(p2, xm);


            if (f1 * fm <= 0) {
                x2 = xm; 
            } else {
                x1 = xm; ;
            }
        }

        return (x1 + x2) / 2; 
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n inner sample points and computing the segment-path between them.
     * assuming x1 < x2.
     * This function should be implemented iteratively (none recursive).
     *
     * @param p                - the polynomial function
     * @param x1               - minimal value of the range
     * @param x2               - maximal value of the range
     * @param numberOfSegments - (A positive integer value (1,2,...).
     * @return the length approximation of the function between f(x1) and f(x2).
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {

        double ans = 0; 

        double dx = (x2 - x1) / numberOfSegments;

        double iznchX = x1;
        double iznchY = Ex1.calPolynom(p, iznchX);

        for (int i = 1; i <= numberOfSegments; i++) {
            double sledX = x1 + i * dx;
            double sledY = Ex1.calPolynom(p, sledX);

            double segment = Math.sqrt(
                    (sledX - iznchX) * (sledX - iznchX) +
                            (sledY - iznchY) * (sledY - iznchY)
            );

            ans += segment;

            iznchX = sledX;
            iznchY = sledY;
        }

        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
     * This function computes an approximation of the area between the polynomial functions within the x-range.
     * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
     *
     * @param p1                - first polynomial function
     * @param p2                - second polynomial function
     * @param x1                - minimal value of the range
     * @param x2                - maximal value of the range
     * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
     * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        double dx = (x2 - x1) / numberOfTrapezoid;

        double iznchX = x1;
        double iznchY = Ex1.calPolynom(p1, iznchX) - Ex1.calPolynom(p2, iznchX);

        for (int i = 1; i <= numberOfTrapezoid; i++) {
            double sledX = x1 + i * dx;
            double sledY = Ex1.calPolynom(p1, sledX) - Ex1.calPolynom(p2, sledX);

            if (iznchY * sledY >= 0) {
                ans += (Math.abs(iznchY) + Math.abs(sledY)) * 0.5 * dx;
            } else {
                double xr = Ex1.sameValue(p1, p2, iznchX, sledX, Ex1.EPS);

                double dxLeft = xr - iznchX;
                double dxRight = sledX - xr;

                double leftArea = 0.5 * Math.abs(iznchY) * dxLeft;
                double rightArea = 0.5 * Math.abs(sledY) * dxRight;

                ans += leftArea + rightArea;
            }

            iznchX = sledX;
            iznchY = sledY;
        }

        return ans;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     */
    public static double[] getPolynomFromString(String p) {
        if (p == null || p.trim().equals("")) {
            return ZERO;
        }
        double[] ans = ZERO;//  -1.0x^2 +3.0x +2.0

        p = p.replace(" ", "");
        p = p.replace("-", "+-");

        if (p.startsWith("+")) {
            p = p.substring(1);
        }
        String[] terms = p.split("\\+");
        int maxDegree = 0;
        for (String term : terms) {
            if (term.equals("")) continue;
            if (term.contains("x^")) {
                int pow = Integer.parseInt(term.substring(term.indexOf("^") + 1));
                if (pow > maxDegree) {
                    maxDegree = pow;
                }
            } else if (term.contains("x")) {
                if (maxDegree < 1) {
                    maxDegree = 1;
                }
            }

        }
        ans = new double[maxDegree + 1];
        for (String term : terms) {
            if (term.equals("")) continue;

            double coef;
            int pow;

            if (term.contains("x^")) {
                String coefStr = term.substring(0, term.indexOf("x"));
                String powStr = term.substring(term.indexOf("^") + 1);

                coef = parseCoefSimple(coefStr);
                pow = Integer.parseInt(powStr);
            } else if (term.contains("x")) {
                String coeffStr = term.substring(0, term.indexOf("x"));

                coef = parseCoefSimple(coeffStr);
                pow = 1;
            } else {
                coef = Double.parseDouble(term);
                pow = 0;
            }
            ans[pow] = coef;
        }
        return ans;
    }

    /**
     * this code helps me distribute the degree correctly
     **/
    private static double parseCoefSimple(String s) {
        if (s.equals("") || s.equals("+")) return 1;
        if (s.equals("-")) return -1;
        return Double.parseDouble(s);
    }


    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     *
     * @param p1
     * @param p2
     * @return
     */
    public static double[] add(double[] p1, double[] p2) {
        double[] ans = ZERO;
        int maxDegree = Math.max(p1.length, p2.length);
        ans = new double[maxDegree];
        for (int i = 0; i < maxDegree; i++) {
            double coef1 = i < p1.length ? p1[i] : 0;
            double coef2 = i < p2.length ? p2[i] : 0;
            ans[i] = coef1 + coef2;
        }
        return ans;
    }

    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     *
     * @param p1
     * @param p2
     * @return
     */
    public static double[] mul(double[] p1, double[] p2) {
        double[] ans = ZERO;
        int degree = (p1.length - 1) + (p2.length - 1);
        ans = new double[degree + 1];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i + j] += p1[i] * p2[j];
            }
        }
        return ans;
    }

    /**
     * This function computes the derivative of the p0 polynomial function.
     *
     * @param po
     * @return
     */
    public static double[] derivative(double[] po) {
        double[] ans = ZERO;
        if (po == null || po.length == 0 || (po.length == 1 && po[0] == 0)){
            return ans;
        }
        if (po.length == 1) {
            return ans;
        }// if degree id zero so a derivative of the number is zero too

        double[] der = new double[po.length - 1];
        for (int i = 1; i < po.length; i++) {
            der[i - 1] = po[i] * i;
        }
		return der;
       
    }
}
