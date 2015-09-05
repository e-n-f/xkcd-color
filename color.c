#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// http://rsb.info.nih.gov/ij/plugins/download/Color_Space_Converter.java

/**
* sRGB to XYZ conversion matrix
*/
double M[3][3] = { {0.4124, 0.3576,  0.1805},
          {0.2126, 0.7152,  0.0722},
          {0.0193, 0.1192,  0.9505} };

/**
* XYZ to sRGB conversion matrix
*/
double Mi[3][3] = { {3.2406, -1.5372, -0.4986},
           {-0.9689,  1.8758,  0.0415},
           { 0.0557, -0.2040,  1.0570} };

double whitePoint[3] = { 95.0429, 100.0, 108.8900 }; /* D65 */

double fpow(double base, double e) {
	return exp(log(base) * e);
}

void RGBtoXYZ(int R, int G, int B, double *x, double *y, double *z) {
	// convert 0..255 into 0..1
	double r = R / 255.0;
	double g = G / 255.0;
	double b = B / 255.0;

	// assume sRGB
	if (r <= 0.04045) {
		r = r / 12.92;
	} else {
		r = fpow(((r + 0.055) / 1.055), 2.4);
	}
	if (g <= 0.04045) {
		g = g / 12.92;
	} else {
		g = fpow(((g + 0.055) / 1.055), 2.4);
	}
	if (b <= 0.04045) {
		b = b / 12.92;
	} else {
		b = fpow(((b + 0.055) / 1.055), 2.4);
	}

	r *= 100.0;
	g *= 100.0;
	b *= 100.0;

	// [X Y Z] = [r g b][M]
	*x = (r * M[0][0]) + (g * M[0][1]) + (b * M[0][2]);
	*y = (r * M[1][0]) + (g * M[1][1]) + (b * M[1][2]);
	*z = (r * M[2][0]) + (g * M[2][1]) + (b * M[2][2]);
}

void XYZtoLAB(double X, double Y, double Z, double *l, double *a, double *b) {
	double x = X / whitePoint[0];
	double y = Y / whitePoint[1];
	double z = Z / whitePoint[2];

	if (x > 0.008856) {
		x = fpow(x, 1.0 / 3.0);
	} else {
		x = (7.787 * x) + (16.0 / 116.0);
	}
	if (y > 0.008856) {
		y = fpow(y, 1.0 / 3.0);
	} else {
		y = (7.787 * y) + (16.0 / 116.0);
	}
	if (z > 0.008856) {
		z = fpow(z, 1.0 / 3.0);
	} else {
		z = (7.787 * z) + (16.0 / 116.0);
	}

	*l = (116.0 * y) - 16.0;
	*a = 500.0 * (x - y);
	*b = 200.0 * (y - z);
}

void LABtoLCH(double L, double A, double B, double *l, double *c, double *h) {
	double C = sqrt(A * A + B * B);
	double H = atan2(B, A);

	*l = L;
	*c = C;
	*h = H;
}

void LCHtoLAB(double L, double C, double H, double *l, double *a, double *b) {
	double A = C * cos(H);
	double B = C * sin(H);

	*l = L;
	*a = A;
	*b = B;
}

void LABtoXYZ(double L, double a, double b, double *xo, double *yo, double *zo) {
	double y = (L + 16.0) / 116.0;
	double y3 = fpow(y, 3.0);
	double x = (a / 500.0) + y;
	double x3 = fpow(x, 3.0);
	double z = y - (b / 200.0);
	double z3 = fpow(z, 3.0);

	if (y3 > 0.008856) {
		y = y3;
	} else {
		y = (y - (16.0 / 116.0)) / 7.787;
	}
	if (x3 > 0.008856) {
		x = x3;
	} else {
		x = (x - (16.0 / 116.0)) / 7.787;
	}
	if (z3 > 0.008856) {
		z = z3;
	} else {
		z = (z - (16.0 / 116.0)) / 7.787;
	}

	*xo = x * whitePoint[0];
	*yo = y * whitePoint[1];
	*zo = z * whitePoint[2];
}

int XYZtoRGB(double X, double Y, double Z, int *R, int *G, int *B) {
	double x = X / 100.0;
	double y = Y / 100.0;
	double z = Z / 100.0;

	// [r g b] = [X Y Z][Mi]
	double r = (x * Mi[0][0]) + (y * Mi[0][1]) + (z * Mi[0][2]);
	double g = (x * Mi[1][0]) + (y * Mi[1][1]) + (z * Mi[1][2]);
	double b = (x * Mi[2][0]) + (y * Mi[2][1]) + (z * Mi[2][2]);

	// assume sRGB
	if (r > 0.0031308) {
		r = ((1.055 * fpow(r, 1.0 / 2.4)) - 0.055);
	} else {
		r = (r * 12.92);
	}
	if (g > 0.0031308) {
		g = ((1.055 * fpow(g, 1.0 / 2.4)) - 0.055);
	} else {
		g = (g * 12.92);
	}
	if (b > 0.0031308) {
		b = ((1.055 * fpow(b, 1.0 / 2.4)) - 0.055);
	} else {
		b = (b * 12.92);
	}

	if (r <= 0 || g <= 0 || b <= 0) {
		return 0;
	}
	if (r > 1 || g > 1 || b > 1) {
		return 0;
	}
	if (isnan(r) || isnan(g) || isnan(b)) {
		return 0;
	}

	*R = r * 255;
	*G = g * 255;
	*B = b * 255;

	return 1;
}

int main(int argc, char **argv) {
	char s[2000];
	while (fgets(s, 2000, stdin)) {
		if (strncmp(s, "INSERT INTO \"answers\"", 16) == 0) {
			char name[2000];
			int r, g, b, id, uid;
			float when;

			if (sscanf(s, "INSERT INTO \"answers\" VALUES(%d,%d,%f,%d,%d,%d,'%[^']'", &id, &uid, &when, &r, &g, &b, name) == 7) {
				double x, y, z;
				RGBtoXYZ(r, g, b, &x, &y, &z);
				double L, A, B;
				XYZtoLAB(x, y, z, &L, &A, &B);
				double c, h;
				LABtoLCH(L, A, B, &L, &c, &h);

				printf("%f %f %f %d %d %d #%02x%0x2%02x %s\n", L, c, h, r, g, b, r, g, b, name);
				// printf("ok: %s", s);
			} else {
				// printf("fail: %s", s);
			}
		}
	}

}
