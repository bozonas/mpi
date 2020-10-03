/* serial.c -- serial trapezoidal rule
 *
 * Calculate definite integral using trapezoidal rule.
 * The function f(x) is hardwired.
 * Input: a, b, n.
 * Output: estimate of integral from a to b of f(x)
 *    using n trapezoids.
 *
 * See Chapter 4, pp. 53 & ff. in PPMPI.
 */

#include <stdio.h>

main() {
    float  integral;   /* Store result in integral   */
    float  a, b;       /* Left and right endpoints   */
    int    n;          /* Number of trapezoids       */
    float  h;          /* Trapezoid base width       */
    float  x;
    int    i;

    float f(float x);  /* Function we're integrating */

    printf("Enter a, b, and n\n");
    scanf("%f %f %d", &a, &b, &n);

    h = (b-a)/n;
    integral = f(a) + f(b);
    x = a;
    for (i = 1; i <= n-1; i++) {
        x = x + h;
        float result = f(x)*2.0;
        if (i%2==1) {
            result*=2.0;
        }
        integral += result;
    }
    integral = integral*h/3.0;

    printf("With n = %d trapezoids, our estimate\n",
        n);
    printf("of the integral from %f to %f = %f\n",
        a, b, integral);
} /* main */


float f(float x) {
    float return_val;
    /* Calculate f(x).  Store calculation in return_val. */

    return_val = x*x;
    return return_val;
} /* f */
