#include <stdio.h>
#include "math.h"
//sqrt function return aquare root of number
double sqrt(double x) {
	double y;
	double low = 0, high = x, b;
	double mid = (low + high) / 2;
	int c = 0;
	if(x < 0)
		//suare root of negative no can't be done , it will gone into imaginary case
		printf("INVALID INPUT");
	else if(x = 0)
		return 0;
		else {
			for (; c <= 50; c++) {
				b = mid * mid;
				if(b == x) {
	        			return mid;//here we get square root
	           			break;
	        		} 
				else {
	        			if(mid * mid > x) {
	        				high = mid;
	                			mid = (low + high) / 2;
	                		} 
					else {
               					low = mid;
               					mid = (low + high) / 2;
            				}
         			}
      			}
		}
}
//pow function returns x^y 
double pow(double x, double y) {
	int i = 1;
	double p;
	if(y == 0)
		return 1;
	if(y > 0) {
		while(i <= y) {
			p = i * x;
			i++;
		}
		return p;
	}
	if(y < 0) {
		int j = 1;
		y = - y;
		while(j <= y) {
			p = j * x;
			j++;
		}
		return (i / p);
	}
}
//fact  function returns factorial value of given number
																																								double fact(unsigned int x) {
	double i = 1;
	double k = 1;
	if (x == 0)
		return 1;
	else {
		for(i = 1; i <= x; i++) {
			k = k * i;
		}
	}
	return k;
}
double sin(double x) {
	//Returns the sine of a given radian angle 
	int n=1;
	double sin = x, t = x;
	while (!(t >= -0.0000000001 && t <= 0.000000001)) {
		t = (-t) * (x * x) / ((2 * n + 1) * (2 * n));
		sin = sin + t;	
		n++;
	}
	return sin;
}
double cos(double x) {
	//Returns the cosine of a given radian angle	
	int n = 1;
	double cos = 1, t = 1;
	while (!(t >= -0.0000000001 && t <= 0.000000001)) {
		t = (-t) * (x * x) / ((2 * n - 1) * (2 * n));
		cos = cos + t;	
		n++;
	}
	return cos;
}
double tan(double x) {
	//Returns the tangent of a given radian angle	
	double tan;
	tan = sin(x) / cos(x);
	return tan;
}
double sinh(double x) {
	//Returns the hyperbolic sine of a given radian angle
	int n = 1;
	double sinh = x, t = x;
	while (!(t >= -0.0000000001 && t <= 0.000000001)) {
		t = (t) * (x * x) / ((2 * n + 1) * (2 * n));
		sinh = sinh + t;	
		n++;
	}
	return sinh;
}
double cosh(double x) {
	//Returns the hyperbolic cosine of a given radian angle
	int n = 1;
	double cosh = 1, t = 1;
	while (!(t >= -1.0000000001 && t <= 1.000000001)) {
		t = (t) * (x * x) / ((2 * n - 1) * (2 * n));
		cosh = cosh + t;	
		n++;
	}
	return cosh;
}
double tanh(double x) {
	//Returns the tangent of a given radian angle
	double tanh;
	tanh =	sinh(x) / cosh(x);
	return tanh;
}
double asin(double x) {
	//Returns the arc (angle) sine of value in radians
	//this is sin inverse function
	double asin = x, t = x;
	int n = 1;
	while (!(t >= -0.0000001 && t <= 0.0000001)) {
		t = ((t) * (x * x) * (2 * n - 1) * (2 * n - 1)) / ((2 * n + 1) * (2 * n));
		asin = asin + t;	
		n++;
	}
	return asin;
}
double acos(double x) {
	//Returns the arc (angle) cosine of value in radians
	//this is cosin inverse function
	
	double acos;
	acos = (PI / 2) - asin(x);
	return acos;
}
double atan(double x) {
	//Returns the arc (angle) sine of value in radians
	//this is sin inverse function
	double atan = x;
	int n = 1;
	double t = x;
	while (!(t >= -0.001 && t <= 0.001)) {
		t = (-t) * (x * x) * (2 * n - 1) / (2 * n + 1);
		atan = atan + t;	
		n++;
	} 
	return atan;
}
double exp(double x) {
	//Returns the value of e raised to the xth power
	double ex = 1;
	double t = 1;
	int n = 1;
	while (!(t >= -0.001 && t <= 0.001)) {
		t = (t * x) / n;
		ex = ex + t;
		n++;
	}	
	return ex;
}

double fabs(double x) {
	//Returns the absolute value of x
	if(x < 0)
		return ( -x);
	else
		return x;
}
double floor(double x) {
	//Returns the largest integer value less than or equal to x
	int y=(int) x;
	int i;
	if (i<x)
		return (double)(++i);
	else 	if(i>x)
			return (double)(i);
		else
			return x;
}		
		
double ceil(double x) {
	//Returns the smallest integer value greater than or equal to x
	int y=(int)x;
	int i;
	if(i<x)
		return (double)(++i);
	else 	if(i>x)
			return (double)(i);
		
		else 
			return x;
}
double log(double x) {
	//Returns base-e logarithm of x
	int n=1,i;
	double log = x - 1;
	int t = 1;
	if (x == 0) 
		printf("value should be greater than zero");
	else {
		while(!(t >= -0.0000000001 && t <= 0.000000001)) {
			t=(-1) * pow((x-1),n) / n;
			log = log +t;
			n = n + 2;
		}
	return log;
}
}
	
double log10(double x) {
	//Returns the base-10 logarithm of x
	double p,log10;
	if (x==0)
		printf("value should be greater than zero");
	
	else { 
		int n=1,i;
		double log = x - 1;
		int t = 1;
		while(!(t >= -0.0000000001 && t <= 0.000000001)) {
			t=(-1) * pow((x-1),n) / n;
			log = log +t;
			n = n + 2;
		}
		p = log;
		log10 = (p / 2.3025850929940);
		return log10;

}
}
double mean(double a[], int n) { // n is total no entered
	double mean,sum=0;
	int i;
	for (i=0;i<n;i++) {
		sum = sum + a[i];
	}
	mean = (float)sum / n;
	return mean;
	}	
double median(double x[], int n) {
    	double tmp;
    	int i, j;
    	// the following two loops sort the array x in ascending order
    	for(i = 1; i < n; i++) {
        for(j=i+1; j < n - 1; j++) {
            if(x[j] < x[i]) {
                // swap elements
                tmp = x[i];
                x[i] = x[j];
                x[j] = tmp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

double mode(double a[], int n) {
	int maxval = 0, maxcount = 0, i, j, count=0;
	for(i = 0; i < n; ++i) {
    		for(j = i + 1; j < n; ++j) {
        		if(a[j] == a[i])
        			++count;
		}
		if (count > maxcount) {
	        	maxcount = count;
         		maxval = a[i];
      		}
   	}

}   	

														






