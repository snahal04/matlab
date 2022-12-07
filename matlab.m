program 1.1 Generation of elementary signals in discrete time

clc; close all;clear all;
%unit impulse sequence
n=-10:1:10;
impulse=[zeros(1,10),ones(1,1),zeros(1,10)];
subplot(2,2,1);stem(n,impulse);
xlabel('discreat time n---->');ylabel("amplitude--->");
title("unit impulse sequence");
axis([-10 10 0 1.2]);
%unit step sequence 
n=-10:1:10;
step=[zeros(1,10),ones(1,11)];
subplot(2,2,2);stem(n,step);
xlabel("discreat time n---->");ylabel("amplitude--->");
title("unit step sequence");
axis([-10 10 0 1.2]);
%unit ramp sequence
n=0:1:10;
ramp=n;
subplot(2,2,3);stem(n,ramp);
xlabel("discreat time n---->");ylabel("amplitude--->");
title("unit ramp sequence");
%unit parabolic sequence
n=0:1:10;
parabola=0.5*(n.^2);
subplot(2,2,4);stem(n,parabola);
xlabel("discreat time n---->");ylabel("amplitude--->");
title("unit parabolic sequence");


program 1.2 Generation of a discrete time exponential sequence

clc;close all;clear all;
n=-10:1:10;
%for 0<a<1
a=0.8;
x1=a.^n;
subplot(2,2,1);stem(n,x1);
title('x1(n)for 0<a<1');
%for a>1
a=1.5;
x2=a.^n;
subplot(2,2,2);stem(n,x2);
title('x2(n)for a>1');
%for -1<a<0
a=-0.8;
x3=a.^n;
subplot(2,2,3);stem(n,x3);
title('x3(n)for -1<a<0');
%for a<-1
a=-1.2;
x4=a.^n;
subplot(2,2,4);stem(n,x4);
title('x4(n)for a<-1');
xlabel('samples n');ylabel('sample amplitude');


program 1.3 Multiplication of discerte-time signal

clc;close all;clear all;
%x1(n)=6*a^n;
n=0:0.1:5;
a=2;
x1=6*(a.^n);
subplot(3,1,1);stem(n,x1);
title('x1(n)');
%x2(n)=2*cos(wn);
f=1.2;
x2=2*cos(2*pi*f*n);
subplot(3,1,2);stem(n,x2);
title('x2(n)');
%multiplication of two sequences
y=x1.*x2;
subplot(3,1,3);stem(n,y);
xlabel('time n');ylabel('amplitude');
title('y(n)');


program 1.4 Even and Odd components of the sequence y(n)=u(n)-u(n-10)


n=-15:1:15;
y1=[zeros(1,15),ones(1,10),zeros(1,6)];
y2=fliplr(y1);
ye=0.5*(y1+y2);
yo=0.5*(y1-y2);
subplot(2,2,1);stem(n,y1);
xlabel('time-->');ylabel('Amplitude-->');
title('y(n)');
subplot(2,2,2);stem(n,y2);
xlabel('time-->');ylabel('Amplitude-->');
title('y(-n)');
subplot(2,2,3);stem(n,ye);
xlabel('time-->');ylabel('Amplitude-->');
title('ye(n)');
subplot(2,2,4);stem(n,yo);
xlabel('time-->');ylabel('Amplitude-->');
title('yo(n)');


program 1.5 Genration of the composite sequence x(n)=u(n+3)+5u(n-15)+4u(n+1)


clc;close all;clear all;
n=-20:1:20;
u=[zeros(1,20),ones(1,21)];
u1=[zeros(1,17),ones(1,24)];
u2=[zeros(1,35),ones(1,6)];u2=5*u2;
u3=[zeros(1,10),ones(1,31)];u3=4*u3;
x=u1+u2+u3;
subplot(4,1,1);stem(n,u1);
title('u(n+3)');
subplot(4,1,2);stem(n,u2);
title('5u(n-15)');
subplot(4,1,3);stem(n,u3);
title('4u(n+10)');
subplot(4,1,4);stem(n,x);
title('x(n)');

program 1.6  Generation of swept frequency sinusoidal signal

clc ; clear all ; close all ;
n = 0 : 100 ; 
a = pi / 2 / 100 ;
b = 0 ;
arg = a*n.*n + b * n ; 
x = cos ( arg ) ; 
stem ( n , x );
xlabel ( ' Discrete time ' ) ;
ylabel ( ' Amplitude ' ) ;
title ( ' Swept - frequency sinusoidal signal ' );

program 1.7 Checking the Time - invariance property 


clc ; clear all;close all ; 
n = 0 : 40 ;
D = 10 ; 
x = 3*cos(2*pi*0.1*n)-2*cos(2*pi*0.4*n) ; 
xd=[zeros(1,D) x] ;
num = [2.2403 2.4908 2.2403] ; 
den=[1 -0.4 0.75]; 
ic = [0 0];
y=filter(num,den,x,ic);
yd = filter ( num , den , xd , ic ) ;
d = y - yd ( 1 + D : 41 + D ) ; 
subplot ( 3,1,1 ) ; stem ( n , y ) ; 
xlabel ( ' Discrete time ' ); ylabel ( "Amplitude " );
title ( ' output y [ n ] ' );
subplot ( 3,1,2 ) ; stem ( n,yd ( 1:41 ) ) ; 
xlabel ( ' Discrete time ' ); ylabel ( ' Amplitude ' );
title ( ' output due to delayed input ' );
subplot ( 3,1,3 ) ; stem ( n , d ) ; 
xlabel ( ' Discrete time ' ) ; ylabel ( ' Amplitude ' );
title ( ' difference signal ' );



program 1.8 Computation of impulse response 

clc ; clear all ; close all ;
N = 40 ;
num = [ 2.2403 2.4908 2.2403 ] ;
den= [ 1 -0.4 0.75 ] ;
y = impz ( num , den , N ) ;
stem ( y ) ; 
xlabel ( ' Discrete time ' );
ylabel ( ' Amplitude ' );
title ( ' Impulse response of the filter ' );


program 1.9 checking the linearity of a system

clc;clear all;close all;
n=0:50;a=2;b=-3;
x1=cos(2*pi*0.1*n);
x2=cos(2*pi*0.4*n);
x=a*x1+b*x2;
num=[2.2403 2.4908 2.2403];
den=[1 -0.4 0.75];
ic=[0 0];
y1=filter(num,den,x1,ic);
y2=filter(num,den,x2,ic);
y=filter(num,den,x,ic);
yt=a*y1+b*y2;
d=y-yt;
subplot(3,1,1);stem(n,y);
xlabel('Discrete time');
ylabel('Amplitude');
title('output due to weighted input');
subplot(3,1,2);stem(n,yt);
xlabel('Discrete time');
ylabel('Amplitude');
title('weighted input');
subplot(3,1,3);stem(n,d);
xlabel('Discrete time');
ylabel('difference signal');

program 1.10 Testing the stability of a system 


clc;clear all; close all;
num = [1 -0.8];
den = [1 1.5 0.9];
N = 200;
h = impz(num,den,N+1);
parsum = 0;
for k = N+1
    parsum = parsum+abs(h(k));
    if abs(h(k))<10^(-6)
        break
    end
end
stem(h);
xlabel('Discrete time');
ylabel('Amplitude');
disp('values =');
disp(abs(h(k)));

program 2.1 Convolution of two sequences


clc;clear all;close all;
x1=[1 2 0 1];
x2=[ 2 2 1 1];
y=conv(x1,x2);
disp('the convol output is');
disp(y)
subplot(3,1,1);
stem(x1);
xlabel('discrete time')
ylabel('amplitude')
title('first input sequence')
subplot(3,1,2);
stem(x2);
xlabel('discrete time')
ylabel('amplitude')
title('second input sequence')
subplot(3,1,3);
stem(y);
xlabel('discrete time')
ylabel('amplitude')
title('convolution output')

program 2.2 Linear convolution via circular convolution 

clc;clear all;close all;
x1=[1 2 3 4 5];
x2=[2 2 0 1 1];
x1e=[x1 zeros(1,length(x2)-1)];
x2e=[x2 zeros(1,length(x1)-1)];
ylin=cconv(x1e,x2e,length(x1e));
disp('linear convolution via circular convolution')
disp('ylin');
y=conv(x1,x2);
disp('Direct convolution');
disp(y);

program 2.3 Linear convolution using DFT

x=[1 2];
h=[2 1];
x1=[x zeros(1,length(h)-1)];
h1=[h zeros(1,length(x)-1)];
X=fft(x1);
H=fft(h1);
y=X.*H;
y1=ifft(y);
disp('the linear conv of the given ')
disp(y1)


program 2.4 Circular convolution using DFT based approach


clc; clear all; close all;
x1=[1 2 0 1];
x2=[2 2 1 1];
d4=[1 1 1 1;1 -j -1 j;1 -1 1 -1;1 j -1 -j];
x11=d4*x1';
x21=d4*x2';
X=x11.*x21;
x=conj((d4)*X/4);
disp('circular convolution by using DFT method')
disp(x') 
x3=cconv(x1,x2,4);
disp('circular convolution by using time domian method')
disp(x3)

program 2.5 Computation of correlation

x1=[1 3 0 4];
y=xcorr(x1,x1);
subplot(2,1,1);stem(x1);
xlabel('Discrete time');
ylabel('Amplitude');
title('Input sequence');
subplot(2,1,2);stem(y);
title('Autocorelation of the input sequence');
xlabel('Discrete time');
ylabel('Amplitude');

program 3.1 Z-transform and inverse Z-transform of given signals


program 3.2 Finding the residues of z^3/((z-0.5)*(z-0.75)*(z-1))


clc;clear all; close all;
syms z;
d=(z-0.5)*(z-0.75)*(z-1);
a1=collect(d);
den=sym2poly(a1);
num=[0 1 0 0];
[num1,den1]=residue(num,den);
fprintf('r1=%4.2f\t',num(1));fprintf('p1=%4.2f\t',den(1));
fprintf('r2=%4.2f\t',num(2));fprintf('p2=%4.2f\t',den(2));
fprintf('r3=%4.2f\t',num(3));fprintf('p3=%4.2f\t',den(3));


program 3.3 Inverse Z-transform by the polynomial division method 
X(Z)=(1+2*z^(-1)+z^(-2))/(1-z^(-1)+0.3561*z^(-2))


clc;clear all;close all;
b=[1 2 1];
a=[1 -1 0.3561];
n=5;
b=[b zeros(1,n-1)];
[h,om]=deconv(b,a);
disp('The terms of inverse z-transforms');
disp(h);


program 3.4 Inverse Z-transform by the polynomial division method 

clc;clear all;close all;
n=5;
n1=[1 -0.22346 1];
n2=[1 -0.4377883 1];
n3=[1 1 0];
d1=[1 -1.433509 0.85811];
d2=[1 -1.293601 0];
d3=[1 -0.612159 0];
b=[n1; n2; n3];
a=[d1; d2; d3];
[b,a]=sos2tf([b a]);
b=[b zeros(1,n-1)];
[x r]=deconv(b,a);
disp('The first five values of inverse z-transform are');
disp(x);
