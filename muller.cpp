//
//  muller.cpp
//
//  Created by WeiKevin on 2015/6/17.
//  Copyright (c) 2015年 WeiKevin. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;

int maxstemp=100; //最高步驟次數
double eps=0.0001; //期望最大誤差值

class complex {
public:
    double x,y;
    complex () {};
    complex (double a,double b) : x(a), y(b) {}
    //標準加法
    complex operator + (const complex&);
    //標準剪法
    complex operator - (const complex&);
    //標準乘法
    complex operator * (const complex&);
    //標準除法
    complex operator / (const complex&);
    //長度專用，假設複數a+bi，則輸出a*a+b*bi
    complex operator *= (const complex&);
};

complex complex::operator + (const complex& param) {
    complex temp;
    temp.x = x + param.x;
    temp.y = y + param.y;
    return temp;
}

complex complex::operator - (const complex& param) {
    complex temp;
    temp.x = x - param.x;
    temp.y = y - param.y;
    return temp;
}

complex complex::operator * (const complex& param) {
    complex temp;
    temp.x = x*param.x - y*param.y;
    temp.y = y*param.x + x*param.y;
    return temp;
}

complex complex::operator / (const complex& param) {
    complex temp;
    temp.x = (x*param.x + y*param.y)/(pow(param.x, 2.0)+pow(param.y, 2.0));
    temp.y = (y*param.x - x*param.y)/(pow(param.x, 2.0)+pow(param.y, 2.0));
    return temp;
}
//複數與常數相乘
complex consttocomplex (double a, complex xcomplex){
    complex answer;
    answer.x = a*xcomplex.x;
    answer.y = a*xcomplex.y;
    return answer;
}

complex complex::operator *= (const complex& param) {
    complex temp;
    temp.x = x * param.x;
    temp.y = y * param.y;
    return temp;
}

//計算虛數的長度
double lengthcomplex(complex root){
    complex temp;
    double length;
    temp=root*=root;
    length=sqrt(temp.x+temp.y);
    return length;
}

//對虛數開根號
complex sqrtcomplex(complex complexsqrt){
    complex answer (0.0,0.0);
    double c=lengthcomplex(complexsqrt);
    answer.x=sqrt((c+complexsqrt.x)/2);
    answer.y=sqrt((c+complexsqrt.y)/2);
    if (complexsqrt.y<0) {
        answer.y=-1*answer.y;
    }
    return answer;
}

//方程式
double func(double x,double* ABC){
    return ABC[0]*pow(x, 2.0)+ABC[1]*pow(x, 1.0)+ABC[2];
}

//方程式輸出虛數解
complex xcomplexfunc(complex xcomplex){
    complex answer (0.0,0.0);
    complex temp[4];
    for (int i=0; i<4; i++) {
        temp[i].x=0.0;
        temp[i].y=0.0;
    }
    temp[3]=consttocomplex(0.0,xcomplex*xcomplex*xcomplex);
    temp[2]=consttocomplex(1.0,xcomplex*xcomplex);
    temp[1]=consttocomplex(0.0,xcomplex);
    temp[0].x=1.0;temp[3].y=0.0;
    
    for (int i=0; i<4; i++) {
        answer=answer+temp[i];
    }
    return answer;
}

//數值運算解 complex
complex answer(double* ABC){
    double square;
    //分子  分母
    complex numerator, denominator;
    //數值解of方程式
    complex answercomplex;
    
    square=ABC[1]*ABC[1]-4*ABC[0]*ABC[2];
    numerator.x=-2*ABC[2];numerator.y=0.0;
    if (square>=0) {
        denominator.x=ABC[1];denominator.y=sqrt(square);
    }else{
        denominator.x=ABC[1];denominator.y=sqrt(-square);
    }
    answercomplex=numerator/denominator;
    return answercomplex;
}

//數值運算解 利用complex ABC
complex answercomplex(complex* ABC){
    complex square;
    //分子  分母
    complex numerator, denominator;
    //數值解of方程式
    complex answercomplex;
    
    square=ABC[1]*ABC[1]-consttocomplex (4.0, ABC[0]*ABC[2]);
    numerator=consttocomplex (-2.0, ABC[2]);
    denominator=ABC[1]+sqrtcomplex(square);
    answercomplex=numerator/denominator;
    return answercomplex;
}

void newABC (complex* ABC,complex* x,complex* fx){
    complex temp1=x[0]-x[2],
            temp2=x[1]-x[2],
            temp3=temp1*temp2*(x[0]-x[1]),
            temp4=fx[1]-fx[2],
            temp5=fx[0]-fx[2];
    ABC[2]=fx[2];
    ABC[1]=((temp1*temp1*temp4)-(temp2*temp2*temp5))/temp3;
    ABC[0]=((temp2*temp5)-(temp1*temp4))/temp3;
}

complex muller(double* x,double* ABC, complex thetrue){
    //誤差計算 complex
    cout<<thetrue.x<<" "<<thetrue.y<<endl;
    complex epscomplex(0.01,0.01);
    double fx[3];
    complex errorcomplex (111.0, 111.0),ABCcomplex[3];
    complex tempanswer;
    for (int i=0; i<3; i++) {
        fx[i]=func(x[i], ABC);
    }
    //由實轉虛
    complex xcomplex[4], fxcomplex[3];
    for (int i=0; i<3; i++) {
        xcomplex[i].x=x[i];xcomplex[i].y=0.0;
        fxcomplex[i].x=fx[i];fxcomplex[i].y=0.0;
    }
    //求新func的ABC係數
    newABC(ABCcomplex,xcomplex,fxcomplex);
    //nextpoint  x3
    tempanswer=answercomplex(ABCcomplex);
    xcomplex[3]=xcomplex[2]+answercomplex(ABCcomplex);
    
    int i=2;
    while (i<maxstemp) {
        cout<<i<<" "<<"曲線趨近的根"<<endl;
        cout<<xcomplex[3].x<<" "<<xcomplex[3].y<<"i"<<endl<<endl;
        if (fabs(xcomplex[3].x-thetrue.x)<eps && fabs(xcomplex[3].y-thetrue.y)<eps) {
            return xcomplex[3];
        }
        for (int i=0; i<3; i++) {
            xcomplex[i]=xcomplex[i+1];
        }
        for (int i=0; i<3; i++) {
            fxcomplex[i]=xcomplexfunc(xcomplex[i]);
        }
        //求新func的ABC係數
        
        newABC(ABCcomplex,xcomplex,fxcomplex);
        
        //nextpoint  x3
        tempanswer=answercomplex(ABCcomplex);
        //cout<<tempanswer.x<<" "<<tempanswer.y<<"i"<<endl;
        xcomplex[3]=xcomplex[2]+answercomplex(ABCcomplex);
        
        i++;
    }
    return errorcomplex;
}

int main(int argc, const char * argv[]) {
    // insert code here...
    complex xrootcomplex;
    complex xrootfunc;
    double ABC[3]={-1,0,1};
    xrootcomplex=answer(ABC);
    cout<<"方程式數值解x^2+1=0的兩個虛數解為下"<<endl;
    cout<<xrootcomplex.x<<"+"<<xrootcomplex.y<<"i"<<" & ";
    cout<<xrootcomplex.x<<-1*xrootcomplex.y<<"i"<<endl;
    cout<<"代回驗證"<<endl;
    xrootfunc=xcomplexfunc(xrootcomplex);
    cout<<xrootfunc.x<<" "<<xrootfunc.y<<"i"<<endl<<endl;    
    cout<<"Müller's Method"<<endl;
    double x[3]={-1.0,0.0,2.0};
    complex answerxroot;
    answerxroot=muller(x, ABC, xrootcomplex);
    cout<<endl<<"答案是："<<endl;
    cout<<answerxroot.x<<" "<<answerxroot.y<<"i"<<endl;
    cout<<answerxroot.x<<" "<<-1*answerxroot.y<<"i"<<endl;
    cout<<"誤差："<<endl;
    cout<<fabs(answerxroot.x-xrootcomplex.x)<<"  ";
    cout <<fabs(answerxroot.y-xrootcomplex.y)<<endl;
    
    return 0;
}

