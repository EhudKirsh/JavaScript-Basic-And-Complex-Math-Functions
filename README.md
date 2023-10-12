List of basic math functions I use in their most efficient JavaScript forms:
```js
'use strict'
/*--------------Array(reduce)Methods--------------*/

/*--------------Statistics--------------*/
// Use Each Of These Directly On An Array, i.e. Func([#,#,#,...])

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

,SumProduct=(Array1,Array2)=>{let Result=0;for(let L=Array1.length;--L>=0;){Result+=Array1[L]*Array2[L]}return Result}

,Correlation=(Array1,Array2)=>{//Pearson Coefficient 'R'
    const n=Array1.length,Sum1=Sum(Array1),Sum2=Sum(Array2)
    return (n*SumProduct(Array1,Array2)-Sum1*Sum2)/Math.sqrt((n*SumProduct(Array1,Array1)-Sum1**2)*(n*SumProduct(Array2,Array2)-Sum2**2))
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

/*--------------Polynomial General Solutions--------------*/
//These find all the non-imaginary solutions for x @y=0

,SolveLinear=(m,c)=>-Number((c/m).toFixed(3))//y=m⋅x+c , m = gradient, c = y-intercept

,SolveQuadratic=(a,b,c)=>{//y=a⋅x²+b⋅x+c, also known as parabolic
    if(a==0){//Must use SolveLinear(b,c) @ a=0, otherwise you'd divide by 0
        return SolveLinear(b,c)
    }else{
        const sqrt_discriminant=Math.sqrt(b**2-4*a*c),a2=2*a
        ,Solutions=[Number((-(sqrt_discriminant+b)/a2).toFixed(3)),Number(((sqrt_discriminant-b)/a2).toFixed(3))].sort((a,b)=>a-b)
        return [Solutions[1],NaN].includes(Solutions[0])?Solutions[0]:Solutions
    }
}
/*--------------Non-Array(reduce)Methods--------------*/
,FracPart=IN=>{ if (String(IN).includes('.') && String(IN).length > 0) { return String(IN).split('.')[1] } else { return '' } }

,RemDiv=(Dividend, Divisor)=>{
    const Sign = Math.sign(Number(Dividend)) * Math.sign(Number(Divisor)); Dividend = String(Dividend).replaceAll('-', ''); Divisor = String(Divisor).replaceAll('-', '')
    const DecFigs = Math.max(FracPart(Dividend).length, FracPart(Divisor).length, 1)/* - 1*/; return Sign * (Dividend % Divisor).toFixed(DecFigs)
}
/*--------------TypeChecks--------------*/
,isNumeric=N=>!isNaN(parseFloat(N))/* && isFinite(n) */
/*--------------Factorial--------------*/
,Gamma=n=>{
    //some magic constants
    const g = 7, // g represents the precision desired, p is the values of p[i] to plug into Lanczos' formula
        p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if(n < 0.5) { return Math.PI / Math.sin(n * Math.PI) / Gamma(1 - n) }
    else { n--; let x = p[0]; for(let i = 1; i < g + 2; i++) { x += p[i] / (n + i); } const t = n + g + 0.5; return Math.sqrt(2 * Math.PI) * Math.pow(t, (n + 0.5)) * Math.exp(-t) * x; }
}
,Factorial=n=>{if (Number.isInteger(n)){if(n>=0){let r = 1; while (n > 0) {r *= n--} return r}else return Infinity} else {return Gamma(n + 1)}}
```
