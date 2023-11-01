List of basic math functions I use in their most efficient JavaScript forms:
```js
'use strict'
/*--------------Array(reduce)Methods--------------*/

/*--------------Statistics--------------*/
// Use Each Of These Directly On An Array, i.e. Func([#,#,...])

const Sum=A=>A.reduce((a,b)=>a+b)

,MeanAverage=A=>Sum(A)/A.length

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

//Sample Variance σ² & Standard Deviation σ
,SampVariance=A=>{const L=A.length,mean=Sum(A)/L;return Sum(A.map(x=>(x-mean)**2))/(L-1)}
,SampSTD=A=>Math.sqrt(SampVariance(A))

//Population Variance S² & Standard Deviation S
,PopVariance=A=>{const L=A.length,mean=Sum(A)/L;return Sum(A.map(x=>(x-mean)**2))/L}
,PopSTD=A=>Math.sqrt(SampVariance(A))

,SumProduct=(A1,A2)=>{let Result=0;for(let L=A1.length;--L>=0;){Result+=A1[L]*A2[L]}return Result}

,Correlation=(A1,A2)=>{//Pearson Coefficient 'R'
    const n=A1.length,Sum1=Sum(A1),Sum2=Sum(A2)
    return (n*SumProduct(A1,A2)-Sum1*Sum2)/Math.sqrt((n*SumProduct(A1,A1)-Sum1**2)*(n*SumProduct(A2,A2)-Sum2**2))
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
,Range=A=>{const a=MinMax3N2(A);return a[1]-a[0]}

,Rank=A=>{const S=A.slice().sort((a,b)=>b-a);return A.map(v=>S.indexOf(v)+1)}
/*--------------Iterative Greatest Common Divisor & Least Common Multiple--------------*/
,gcd=(a,b)=>{a<b&&([a,b]=[b,a]);while(b!=0){[a,b]=[b,a%b]}return a}
,GCD=A=>A.reduce(gcd)

,lcm=(a,b)=>a*b/gcd(a,b)
,LCM=A=>A.reduce(lcm)

/*--------------Factorial !--------------*/
,Gamma=n=>{//The use of this 'Gamma' is to ensure Factorial can work on fractions, not just integers
    const g = 7, // g represents the precision desired, p is the values of p[i] to plug into Lanczos' formula
    //some magic constants
        p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if(n < 0.5) { return Math.PI / Math.sin(n * Math.PI) / Gamma(1 - n) }
    else { n--; let x = p[0]; for(let i = 1; i < g + 2; i++) { x += p[i] / (n + i); } const t = n + g + 0.5; return Math.sqrt(2 * Math.PI) * Math.pow(t, (n + 0.5)) * Math.exp(-t) * x; }
}
,Factorial=n=>{if(Number.isInteger(n)){if(n>=0){let r = 1; while (n > 0) {r *= n--} return r}else return Infinity} else {return Gamma(n + 1)}}

/*--------------Polynomial General Solutions--------------*/
//These find all the real and complex roots for x @ y=0 and real/non-complex coefficients a,b,c,d & e

,SolveLinear=(m,c)=>//y=mx+c , m = gradient, c = y-intercep
    m==0?console.log('No Gradient? No Root!'):console.log('Root (y=0): x = '+Number((-c/m).toFixed(3)))

,SolveQuadratic=(a,b,c)=>{//y=ax²+bx+c, also known as parabolic
    if(a==0){//Must use SolveLinear(b,c) @ a=0, otherwise you'd divide by 0
        return SolveLinear(b,c)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_2
        const discriminant=b**2-4*a*c,a2=2*a,ba2=Number((-b/a2).toFixed(3))
        let Nature='Minima';a<0&&(Nature='Maxima') // δ²y/δx² = 2a
        if(discriminant==0){//1 'repeated' real root which is also the stationary point
            console.log('Root (y=0) & '+Nature+': x = '+ba2)
        }else{//Either 2 real roots OR 2 complex roots
            const sqrt_discriminantOver_2a=Number((Math.sqrt(Math.abs(discriminant))/a2).toFixed(3))
            ,x=Number((-b/(2*a)).toFixed(3)) // δy/δx = 0 = 2ax+b ∴ x=-b/(2a)
            if(discriminant>0){//2 real roots
                const Roots=[Number((ba2+sqrt_discriminantOver_2a).toFixed(3)),Number((ba2-sqrt_discriminantOver_2a).toFixed(3))].sort((a,b)=>a-b)
                console.log('Roots (y=0): x₀ = '+Roots[0]+' , x₁ = '+Roots[1]+'\n'+Nature+': x = '+x+' , y = '+Number((a*x**2+b*x+c).toFixed(3)))
            }else{//2 complex roots
                let Roots=ba2+' ± '+sqrt_discriminantOver_2a+' i';sqrt_discriminantOver_2a==1&&(Roots=ba2+' ± i')
                console.log('Roots (y=0): x = '+Roots+'\n'+Nature+': x = '+x+' , y = '+Number((a*x**2+b*x+c).toFixed(3)))
            }
        }
    }
}
,SolveCubicDiscriminant=(a,b,c,d)=>{//y=a⋅x³+b⋅x²+c⋅x+d
    if(a==0){
        return SolveQuadratic(b,c,d)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_3
        const discriminant=18*a*b*c*d-4*d*b**3+b**2*c**2-4*a*c**3-27*a**2*d**2
        //Depressed cubic, source: https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
        if(discriminant==0){
            return '2 real roots, 1 of which is a stationary point'
        }else if(discriminant>0){
            return '3 real distinct roots'
        }else{
            return '1 real root and 2 complex roots'
        }
    }
}
,SolveQuarticDiscriminant=(a,b,c,d,e)=>{//y=a⋅x⁴+b⋅x³+c⋅x²+d⋅x+e
    if(a==0){
        return SolveCubicDiscriminant(b,c,d,e)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_4
        const discriminant=256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*c**2*e**2+144*a**2*c*d**2*e-27*a**2*d**4+144*a*b**2*c*e**2-6*a*b**2*d**2*e-80*a*b*c**2*d*e+18*a*b*c*d**3+16*a*c**4*e-4*a*c**3*d**2-27*b**4*e**2+18*b**3*c*d*e-4*b**3*d**3-4*b**2*c**3*e+b**2*c**2*d**2
        if(discriminant==0){
            return 'at least 2 roots are equal'
        }else if(discriminant>0){
            return 'Roots are either all real or all complex'
        }else{
            return '2 real roots and 2 complex roots'
        }
    }
}
```
