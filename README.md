List of basic and complex math functions I use in their most efficient JavaScript forms:
```js
'use strict'
/*--------------Statistics--------------*/
// Use Each Of These Directly On An Array, i.e. Func([#,#,...])

const Sum=A=>A.reduce((a,b)=>a+b) // Σ

,SumProduct=(A1,A2)=>{let Π=0;for(let L=A1.length;--L>=0;){Π+=A1[L]*A2[L]}return Π} // ∏

,MeanAverage=A=>Sum(A)/A.length // µ

,MedianAverage=A=>{A.sort((a,b)=>a-b);const L2=A.length/2;return L2%1==0?/*If Even*/(A[L2]+A[L2-1])/2:/*If Odd*/A[L2-0.5]}

,ArrayOccurancesEachItem=A=>{
    const counter={};for(let L=A.length;--L>=0;){const element=A[L];counter[element]?counter[element]+=1:counter[element]=1}return counter
}
,ModeAverage=A=>{//Array items with higest occurances, not just numbers but also strings
    Array.isArray(A)&&(A=ArrayOccurancesEachItem(A))//Input can be array [] or object {}
    const HighestOccurance=Object.values(A).reduce((a,b)=>Math.max(a,b))
    ,keys=Object.keys(A)
    let Mode=filteredKeys=keys.filter(key=>A[key]==HighestOccurance)
    if(Mode.length<=1)Mode=Mode[0]
    console.log('Mode: ',Mode,' Occurance: '+HighestOccurance)
    return [Mode,HighestOccurance]
}

,Minimum=A=>A.reduce((a,b)=>Math.min(a,b))
,Maximum=A=>A.reduce((a,b)=>Math.max(a,b))
,MinMax3N2=A=>{//This is a 25% faster/most efficient way to get both the Minimum and Maximum of an array at the same time
    let L=A.length,Min,Max
    if(L&1){Min=Max=Number(A[--L])}else{const a=Number(A[--L]),b=Number(A[--L]);a>b?(Max=a,Min=b):(Max=b,Min=a)}
    for(let Big,Small;L>0;){
        const a=Number(A[--L]),b=Number(A[--L])
        a>b?(Big=a,Small=b):(Big=b,Small=a)
        Big>Max&&(Max=Big);Small<Min&&(Min=Small)
    }
    return [Min,Max]
}
,Range=A=>{const a=MinMax3N2(A);return a[1]-a[0]} // Max-Min

,Rank=A=>{const S=A.slice().sort((a,b)=>b-a);return A.map(v=>S.indexOf(v)+1)}
// Purpose: rank values from biggest to smallest, e.g. Rank([0, -1.4, 7.6]) => [2, 3, 1]

,FractionRank=A=>{// Adds 0.5*(occurances-1) per each tied rank
    const R=Rank(A),L=R.length,rankCounts=new Map()
    for (let i=-1;++i<L;){
        const rank=R[i]
        rankCounts.has(rank)||rankCounts.set(rank,0)
        rankCounts.set(rank,rankCounts.get(rank)+1)
    }
    for (let i=-1;++i<L;)
        R[i]+=(rankCounts.get(R[i])-1)/2
    return R
}//e.g. FractionRank([5,62,5,0,4,-3.4,5,62]) => [4, 1.5, 4, 7, 6, 8, 4, 1.5]

,SpearmanCorrelation=(A1,A2)=>{// Coefficient 'ρ': 1 ≥ ρ ≥ -1
    const n=A1.length,Rank1=FractionRank(A1),Rank2=FractionRank(A2);let Σd2=0
    for(let L=n;--L>=0;)Σd2+=(Rank1[L]-Rank2[L])**2
    return Number((1-(6*Σd2/(n*(n**2-1)))).toFixed(3))
}
,PearsonCorrelation=(A1,A2)=>{// Coefficient 'R': 1 ≥ R ≥ -1
    const n=A1.length
    if(n<2)return -Infinity // to prevent error in console when used in some apps that may not accept NaN
    const Sum1=Sum(A1),Sum2=Sum(A2)
    return (n*SumProduct(A1,A2)-Sum1*Sum2)/Math.sqrt((n*SumProduct(A1,A1)-Sum1**2)*(n*SumProduct(A2,A2)-Sum2**2))
}

//Population Variance σ² & Standard Deviation σ
,VAR_P=A=>{const L=A.length,µ=Sum(A)/L;return Sum(A.map(x=>(x-µ)**2))/L}
,STDEV_P=A=>Math.sqrt(VAR_P(A))
,COVAR_P=(A1,A2)=>{
    const L=A1.length,µ1=Sum(A1)/L,µ2=Sum(A2)/L;let T=0
    for(let i=L;--i>=0;)T+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((T/L).toFixed(3))
}
//Sample Variance S² & Standard Deviation S
,VAR_S=A=>{const L=A.length,µ=Sum(A)/L;return Sum(A.map(x=>(x-µ)**2))/(L-1)}
,STDEV_S=A=>Math.sqrt(VAR_S(A))
,COVAR_S=(A1,A2)=>{
    const L=A1.length,µ1=Sum(A1)/L,µ2=Sum(A2)/L;let T=0
    for(let i=L;--i>=0;)T+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((T/(L-1)).toFixed(3))
}

,BestFitLine=(X,Y)=>{//Least Square Method: inputs are arrays of x & y coordinates which must be of equal lengths
    let L=X.length,N=D=0
    if(L!=Y.length)return 'Please enter arrays of x and y coordinates of equal lengths'

    const µ_Y=Y.reduce((a,b)=>a+b)/L,µ_X=X.reduce((a,b)=>a+b)/L // Mean averages of each array

    while(--L>=0){
        const x=X[L]-µ_X//y=Y[L]-µ_Y
        N+=x*(Y[L]-µ_Y);D+=x**2
    }
    const m=Number((N/D).toFixed(3))
    return [m,Number((µ_Y-m*µ_X).toFixed(3))]
    // Outputs [m,c] (m = gradient & c = y-intercep of y=mx+c)
}
,LinearExtrapolation=(x,x0,x1,y0,y1)=>Number((y0+(x-x0)*(y1-y0)/(x1-x0)).toFixed(3)) // Outputs y, also works for extrapolation
// Use LinearExtrapolation when you have only 2 other coordinates, otherwise use BestFitLine and y=mx+c

,TrapeziumRule=(X_Interval,Y_Values_Array)=>Number(((2*Sum(Y_Values_Array)-Y_Values_Array[0]-Y_Values_Array.at(-1))*X_Interval/2).toFixed(3))

/*--------------Iterative Greatest Common Divisor & Least Common Multiple--------------*/
,gcd=(a,b)=>{a<b&&([a,b]=[b,a]);while(b!=0){[a,b]=[b,a%b]}return a}
,GCD=A=>A.reduce(gcd)

,lcm=(a,b)=>a*b/gcd(a,b)
,LCM=A=>A.reduce(lcm)

/*--------------Logarithms & Exponentials--------------*/
,BaseLog=(Base,Num)=>Number((Math.log(Num)/Math.log(Base)).toFixed(3))

/*--------------Factorial !--------------*/
,Gamma=n=>{//The use of this 'Gamma' is to ensure Factorial can work on fractions, not just integers
    const g = 7, // g represents the precision desired, p is the values of p[i] to plug into Lanczos' formula
    //some magic constants
        p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if(n < 0.5) { return Math.PI / Math.sin(n * Math.PI) / Gamma(1 - n) }
    else { n--; let x = p[0]; for(let i = 1; i < g + 2; i++) { x += p[i] / (n + i); } const t = n + g + 0.5; return Math.sqrt(2 * Math.PI) * Math.pow(t, (n + 0.5)) * Math.exp(-t) * x; }
}
,Factorial=n=>{if(Number.isInteger(n)){if(n>=0){let r = 1; while (n > 0) {r *= n--} return r}else return Infinity} else {return Gamma(n + 1)}}

/*--------------Imaginary & Complex Numbers--------------*/
,ComplexSQRT=(a,b)=>{// Square Root √(a+bi), for real a & b
    if(b==0){
        return Math.sqrt(a)
    }else{
        const HYPOT_ab=Math.hypot(a,b)
        return [Math.sqrt((HYPOT_ab+a)/2),Math.sign(b)*Math.sqrt((HYPOT_ab-a)/2)] //±[](0) is real, ±[](1) is imaginary
    }
}
,ComplexCBRT=(a,b)=>{// Cube Root ₃√(a+bi), for real a & b
    if(b==0){
        return Math.cbrt(a)
    }else{
        const θ=Math.atan2(b,a)/3,r=Math.cbrt(Math.hypot(a,b))
        return [r*Math.cos(θ),r*Math.sin(θ)] //[](0) is real, [](1) is imaginary
    }
}

/*---Polynomials: General Solutions To Find Roots, Stationary Points & Their Natures---*/
//These find all the real and complex roots for x @ y=0 and stationary points and their natures for real/non-complex coefficients a, b, c, d, e & f

,SolveLinear=(m,c)=>//y=mx+c , m = gradient, c = y-intercep
    m==0?console.log('No Gradient? No Root!'):console.log('Root (y=0): x = '+Number((-c/m).toFixed(3)))

,SolveQuadratic=(a,b,c)=>{//y=ax²+bx+c, also known as parabolic
    if(a==0){//Must use SolveLinear(b,c) @ a=0, otherwise you'd divide by 0
        return SolveLinear(b,c)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_2
        const Discriminant=b**2-4*a*c,a2=2*a,aNegative=a<0,x=Number((-b/a2).toFixed(3))// δy/δx = 0 = 2ax+b ∴ x=-b/(2a)
        let Nature='Minima';aNegative&&(Nature='Maxima') // δ²y/δx² = 2a
        if(Discriminant==0){//1 'repeated' real root which is also the stationary point
            console.log('Root (y=0) & '+Nature+': x = '+x)
        }else{//Either 2 real roots OR 2 complex roots
            const sqrt_discriminantOver_2a=Number((Math.sqrt(Math.abs(Discriminant))/a2).toFixed(3))
            ,y=Number((a*x**2+b*x+c).toFixed(3))
            if(Discriminant>0){//2 real roots
                let x0=Number((x-sqrt_discriminantOver_2a).toFixed(3)),x1=Number((x+sqrt_discriminantOver_2a).toFixed(3)),flip
                if(aNegative){flip=x0;x0=x1;x1=flip} //sort so x₀<x₁
                console.log('• Roots (y=0): x₀ = '+x0+' , x₁ = '+x1+'\n• '+Nature+': x = '+x+' , y = '+y)
            }else{//2 complex roots
                let Roots=x+' ± '+Math.abs(sqrt_discriminantOver_2a)+' i';Math.abs(sqrt_discriminantOver_2a)==1&&(Roots=x+' ± i')
                console.log('• Roots (y=0): x = '+Roots+'\n• '+Nature+': x = '+x+' , y = '+y)
            }
        }
    }
}
,CubicDiscriminant=(a,b,c,d)=>{//y=ax³+bx²+cx+d
    if(a==0){
        return SolveQuadratic(b,c,d)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_3
        const Discriminant=18*a*b*c*d-4*d*b**3+b**2*c**2-4*a*c**3-27*a**2*d**2
        if(Discriminant==0){
            return 'At least 1 stationary point is also a root'
        }else if(Discriminant>0){//∴ Quadratic_Discriminant 4*(b**2-3*a*c)>0
            return '3 real distinct roots'
        }else{
            return "Only 1 real root which isn't a stationary point"
        }
    }
}
,QuarticDiscriminant=(a,b,c,d,e)=>{//y=ax⁴+bx³+cx²+dx+e
    if(a==0){
        return CubicDiscriminant(b,c,d,e)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_4
        const Discriminant=256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*c**2*e**2+144*a**2*c*d**2*e-27*a**2*d**4+144*a*b**2*c*e**2-6*a*b**2*d**2*e-80*a*b*c**2*d*e+18*a*b*c*d**3+16*a*c**4*e-4*a*c**3*d**2-27*b**4*e**2+18*b**3*c*d*e-4*b**3*d**3-4*b**2*c**3*e+b**2*c**2*d**2
        if(Discriminant==0){
            return 'At least 1 stationary point is also a root'
        }else if(Discriminant>0){
            return 'Roots are either all real or all complex'
        }else{
            return '2 real roots and 2 complex roots'
        }
    }
}
```
