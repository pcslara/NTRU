#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


typedef struct __poly {
    short * p;              // p_i coeficientes
    int size;               // N
} poly;


void print_latex( poly, char *);
int mod( int , int );
void print_latex( poly p, char * name ) {
    int i;
    int printSignal = 0;
    printf("%s(x)=", name );
    
    for( i = p.size-1; i >= 0; i-- ) {
        if( p.p[i] == 0 )
            continue;    
        if( i == 0 ) { 
            if( p.p[i] > 0 ) {
                if( printSignal )
                    printf("+");
                else
                    printSignal = 1;
            }
            printf("%d", p.p[i]);
        } else if( i == 1 ) {
            if( p.p[i] > 0 )
                printf("+");
            if( p.p[i] != 1 )
                printf("%dx", p.p[i]);
            else
                printf("x");
        } else {
            if( p.p[i] > 0 ) {
                if( printSignal )
                    printf("+");
                else
                    printSignal = 1;
            }
            if( p.p[i] != 1 )
                printf("%dx^{%d}", p.p[i], i);
            else
                printf("x^{%d}", i);
        }
    }
    printf("\n");

}
poly new_poly( int N ) {
    int i;
    poly pol;
    pol.p = (short *) malloc( sizeof(short) * N );
    pol.size = N;
    for(i = 0; i < N; i++ ) {
        pol.p[i] = 0;
    }
    return pol;
}

void free_poly( poly p ) {
    free( p.p );
}
/**
 * return 1, if p(x)=1
 * return 0, otherwise
 */
short isone( poly p) {
    int i;
    for( i = p.size - 1; i > 0; i-- )
        if( p.p[i] != 0 )
            return 0;
    if( p.p[0] == 1 )
        return 1;
    else
        return 0;
}
/**
 * return 1, if p(x)=0
 * return 0, otherwise
 */
short iszero( poly p ) {
    int i;
    for( i = 0; i  < p.size; i++ )
        if( p.p[i] != 0 )
            return 0;
    return 1;
}

void print( poly p ) {
    int i;
    printf("[");
    for( i = 0; i < p.size; i++ ) {
        printf("%d ", p.p[i] );
    }
    printf("]\n");
}

void zero( poly p ) {
    int i;
    for( i = 0; i < p.size; i++ ) 
        p.p[i] = 0;
}

void copy( poly r, poly a ) {
    int i;
    int size = r.size;
    if( r.size > a.size  ) {
        size = a.size;
    }

    for( i = 0; i < size; i++ )
        r.p[i] = a.p[i];
}

void copy_mod( poly r, poly a, int p) {
    int i;
    int size = r.size;
    if( r.size > a.size  ) {
        size = a.size;
    }

    for( i = 0; i < size; i++ )
        r.p[i] = mod( a.p[i], p );
}

int mod(int x, int N){
    return ((x % N) + N) % N;
}
/**
 * return 1, if a(x)==b(x)
 * return 0, otherwise
 */
short iseq( poly a, poly b) {
    int i;
    if( a.size != b.size  ) {
        printf("DIFF SIZE EQ\n");
        return 0;
    }
    for( i = 0; i  < a.size; i-- )
        if( a.p[i] != b.p[i] )
            return 0;
    return 1;
}

void exchange(poly a, poly b ) {
    int i;
    short tmp;
    if( a.size != b.size  ) {
        printf("DIFF SIZE EXCHANGE\n");
        return;
    }
    for( i = 0; i < b.size; i++ ) {
        tmp = a.p[i];
        a.p[i] = b.p[i];
        b.p[i] = tmp;
    }   
}
/**
 *  f = f << 1 
 */   
 
void shiftRotateL(poly f) {
    int i;
    short fn = f.p[f.size - 1];
    for (int i=f.size-1; i>0; i--)
        f.p[i] = f.p[i-1];
    f.p[0] = fn;
}

/**
 *  f = f >> 1 
 */  
 
void shiftRotateR(poly f) {
    int i;
    short f0 = f.p[0];
    for (int i=0; i < f.size - 1; i++)
        f.p[i] = f.p[i+1];
    f.p[f.size - 1] = f0;
}

 
int deg( poly r ) {
    int i;
    for( i = r.size - 1; i != 0; i-- )
        if( r.p[i] != 0 )
            return i;
    return 0;
}

// mpz_invert

int invmod(int a, int b) {
    if( a == 0 )
        return 0;
	int b0 = b, t, q;
	int x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
        if( b == 0 )
            return 0;
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}
/**
 * r = a + b
 */
 
void add_poly( poly r, poly a, poly b, int n ) {
    int i;
    if( r.size != a.size || r.size != b.size ) {
        printf("DIFF SIZE ADD\n");
        return;
    }
    for( i = 0; i < a.size; i++ )
        r.p[i] = mod(a.p[i] + b.p[i], n );
}
/**
 * r = a - b
 */
void sub_poly( poly r, poly a, poly b, int n ) {
    int i;
    if( r.size != a.size || r.size != b.size ) {
        printf("DIFF SIZE SUB\n");
        return;
    }
    for( i = 0; i < a.size; i++ )
        r.p[i] = mod(a.p[i] - b.p[i], n);
}


void convolution( poly r, poly a, poly b, int n ) {
    int k, j;
    
    if( r.size != a.size || r.size != b.size ) {
        printf("DIFF SIZE CONV\n");
        return;
    }
    for( k = 0; k < r.size; k++ )
        r.p[k] = 0;

    for( k = 0; k < r.size; k++ ){
        for( j = 0; j < r.size; j++ )
            r.p[(k+j) % r.size] += mod( a.p[k] * b.p[j], n ); 
    }

    for( k = 0; k < r.size; k++ ){
        r.p[k] = mod( r.p[k], n );
    }
}

int inverse_prime( poly r, poly a, int p ) {
    int N = a.size;

    // Initialization:
    // k=0, b(X) = 1, c(X) = 0, f(X)=a(X), g(X)=X^N-1
    int k = 0;
    poly b = new_poly(N+1); b.p[0] = 1;
    poly c = new_poly(N+1); 
    poly f = new_poly(N+1); copy_mod( f, a, p );
    poly g = new_poly(N+1); g.p[N] = 1; g.p[0] = (p-1);

    while( 1 ) {
        
        
        while ( f.p[0] == 0  && deg( f ) > 0 ) {
            shiftRotateR(f);
            shiftRotateL(c);
            k++;
       }
        if ( deg( f ) == 0) {
            int f0Inv = invmod( f.p[0], p );
            if (f0Inv == 0)
                return 0;
            int shift = mod( N-k, N );
            for (int i=0; i<N; i++)
                r.p[(i+shift) % N] = mod(f0Inv * b.p[i], p);
            return 1;
        }

        if( deg( f ) < deg( g ) ) {
            exchange( f, g ); exchange( b, c );
        }

        int g0Inv = invmod( g.p[0], p );
        if( g0Inv == 0 )
            return 0;
        
        short u = mod( f.p[0] * g0Inv, p );

        for( int i = 0; i < f.size; i++ ) {
            f.p[i] = mod( f.p[i] - u * g.p[i], p );    
        }

        for( int i = 0; i < b.size; i++ ) {
            b.p[i] = mod( b.p[i] - u * c.p[i], p );    
        }
    }
    free_poly( b );
    free_poly( c );
    free_poly( f );
    free_poly( g );
}

int ipow( int b, int e ) {
    int i;
    int ret = 1;
    for(i = 0; i < e; i++ ) {
        ret *= b;
    }
    return ret;
}

void multiply_constant( poly f, short constant, int n ) {
    int i;
    for(i = 0; i < f.size; i++ ) {
        f.p[i] = mod( constant * f.p[i], n ); 
    }
}


void random_poly_L( poly f, int d1, int d2 ) {
    int i;
    
    for( i = 0; i < f.size; i++ )
        f.p[i] = 0;
    for( i = 0; i < d1; i++ ) {
        int idx;
        do {
            idx = rand() % f.size; 
        } while( f.p[idx] != 0 );

        f.p[idx] = 1;
    }


    for( i = 0; i < d2; i++ ) {
        int idx;
        do {
            idx = rand() % f.size; 
        } while( f.p[idx] != 0 );

        f.p[idx] = -1;
    }

}

int inverse_power_prime( poly b, poly a, int p, int r ) {
    int q = p;
   
    int  pr = ipow( p, r );
    int i;
    
    if( !inverse_prime( b, a, p ) ) {
        return 0;
    }
    
    poly t = new_poly( a.size );
    poly b2 = new_poly(a.size );
    
    while( q < pr ) {
        q = q * q;
        convolution(t, a, b, q);
        for(i = 0; i < a.size; i++ ) {
            t.p[i] = mod( -t.p[i], q );
        }
        t.p[0] = mod( t.p[0] + 2, q );
        if( q > pr ) 
            convolution( b2, b, t, pr );
        else
            convolution( b2, b, t, q );
        copy( b, b2 );    
        
    }
    free_poly( t );
    free_poly( b2 );
    return 1;
}

void logprime( int * base, int * exponent, int q ) {
    int primes [] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 31 };
    int p = 0;
    int i;

    for( i = 0; i < 10; i++ ) {
        if( q % primes[i] == 0 ) {
            p = primes[i];
            break;
        }
    }

    while( q != 1 ) {
        q = q / p;
        i++; 
    }
    *base = p;
    *exponent = i;
}

void center( poly f, int q ) {
    int i;
    for(i = 0; i < f.size; i++ ) {
        if( f.p[i] > q/2 )
            f.p[i] =  mod( f.p[i], q ) - q;
        else
            f.p[i] =  mod( f.p[i], q );
    }
}

void keygen_2005( poly h, poly f, poly fq,  poly F, poly g, int N, int p, int q, int df, int dg, int dF ) {
    
    random_poly_L( g, dg, 0 );
    while( 1 ) {
        random_poly_L( F, dF, 0 );
        copy( f, F );
        for( int i = 0; i < N; i++ ) f.p[i] *= p;
        f.p[0] = f.p[0] + 1;
        if( inverse_prime( fq, f, q  ) == 1  ) {
            break;
        }
    }
    convolution( h, g, fq, q );
    multiply_constant( h, p, q );
}

poly encrypt_2005( poly m, poly h, int p, int q, int dr ) {
    poly c = new_poly( m.size );
    poly t = new_poly( m.size );
    poly R = new_poly( m.size );
    random_poly_L( R, dr, dr );
    copy( c, R );                          // c =  R
    convolution( t, c, h, q );             // t =  R * h mod q
    add_poly( c, t, m, q );                // c =  R * h + m mod q
    return c;
}

poly decrypt_2005( poly c, poly f, int p, int q ) {
    poly a = new_poly( c.size );
    poly m = new_poly( c.size );
    convolution( a, f, c, q );
    center( a, q  );
    convolution( m, f, a, p );
    center( m, p  );
    return m;
}


int main(int argc, char ** argv ) {
    srand( time(NULL) );
    int N = 48;
    int p = 2;
    int q = 173;

    int df = 12;
    int dg = df;
    int dr = df; 
    int dF = df;

    poly h  = new_poly( N );
    poly f  = new_poly( N );
    poly fp = new_poly( N );
    poly fq = new_poly( N );
    poly m = new_poly( N );
    poly F = new_poly( N );
    poly g = new_poly( N );

    m.p[0] = 1;
    m.p[3] = 1;
    m.p[7] = 1;
    
    print( m );
    keygen_2005( h, f, fq, F, g, N, p, q, df, dg, dF );
    
    poly c = encrypt_2005(m, h, p, q, dr );
    print( c );    
    poly m2 = decrypt_2005( c, f, p, q );
    print( m2 );
           
    return 0;
}

