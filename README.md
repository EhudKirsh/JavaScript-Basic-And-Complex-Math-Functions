List of basic and complex math functions I use in their most efficient JavaScript forms:
```js
'use strict'
/*--------------Statistics--------------*/
// Use Each Of These Directly On An Array, i.e. Func([#,#,...])

const Sum=A=>{let L=A.length,S=0;do{S+=A[--L]}while(L>0);return S} // Σ

,Product=A=>{let L=A.length,P=A[--L];while(--L>=0){const a=A[L];if(a==0)return 0;P*=a}return P} // ∏
// Efficiently returns 0 as soon as any is detected. This can also be modified to detect Infinity, -Infinity and NaN

,SumProduct=(A1,A2)=>{let L=A1.length,R=A1[--L]*A2[L];while(--L>=0)R+=A1[L]*A2[L];return R}

,MeanAverage=A=>{const L=A.length;let l=L,Σ=A[--l];while(--l>=0)Σ+=A[l];return Σ/L} // µ

,MedianAverage=A=>{A.sort((a,b)=>a-b);const L2=A.length/2;return L2%1?A[L2-0.5]:(A[L2]+A[L2-1])/2}

,CountOccurances=A=>{const C={};for(let L=A.length;--L>=0;){const E=A[L];C[E]?C[E]+=1:C[E]=1}return C}
,ModeAverage=A=>{//Array items with higest occurances, not just numbers but also strings
    Array.isArray(A)&&(A=CountOccurances(A))//Input can be array [] or object {}
    const HighestOccurance=Maximum(Object.values(A))
    ,keys=Object.keys(A)
    let Mode=filteredKeys=keys.filter(key=>A[key]==HighestOccurance)
    if(Mode.length<2)Mode=Mode[0]
    console.log('Mode: ',Mode,' Occurance: '+HighestOccurance)
    return [Mode,HighestOccurance]
}

,Minimum=A=>{let L=A.length,M=A[--L];while(--L>=0){const N=A[L];N<M&&(M=N)};return M}
,Maximum=A=>{let L=A.length,M=A[--L];while(--L>=0){const N=A[L];N>M&&(M=N)};return M}
,MinMax3N2=A=>{//This is a slightly faster algorithm to get both the Minimum and Maximum of an array at the same time
    let L=A.length,Min,Max
    if(L&1){Min=Max=A[--L]}else{const a=A[--L],b=A[--L];a>b?(Max=a,Min=b):(Max=b,Min=a)}
    for(let Big,Small;L>0;){
        const a=A[--L],b=A[--L]
        a>b?(Big=a,Small=b):(Big=b,Small=a)
        Big>Max&&(Max=Big);Small<Min&&(Min=Small)
    }
    return [Min,Max]
}
,Range=A=>{const a=MinMax3N2(A);return a[1]-a[0]} // Max-Min

,Rank=A=>{const S=A.slice().sort((a,b)=>b-a);return A.map(v=>S.indexOf(v)+1)}
// Purpose: rank values from biggest to smallest, e.g. Rank([0, -1.4, 7.6]) => [2, 3, 1]

,FractionRank=A=>{// Adds 0.5*(occurances-1) per each tied rank
    const R=Rank(A),L=R.length,C=new Map()
    for(let i=-1;++i<L;){const r=R[i];C.has(r)||C.set(r,0);C.set(r,C.get(r)+1)}
    for(let i=-1;++i<L;)R[i]+=(C.get(R[i])-1)/2
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
    let i=n,Sum1=A1[--i],Sum2=A2[i],SumProduct11=Sum1**2,SumProduct22=Sum2**2,SumProduct12=Sum1*Sum2
    while(--i>=0){
        const a1=A1[i],a2=A2[i]
        Sum1+=a1;Sum2+=a2;SumProduct11+=a1**2;SumProduct22+=a2**2;SumProduct12+=a1*a2
    }
    return (n*SumProduct12-Sum1*Sum2)/Math.sqrt((n*SumProduct11-Sum1**2)*(n*SumProduct22-Sum2**2))
}

//Population Variance σ², Standard Deviation σ & Covariance
,VAR_P=A=>{
    const L=A.length
    let i=L,µ=A[--i];while(--i>=0)µ+=A[i];µ/=L
    i=L;let V=(A[--i]-µ)**2;while(--i>=0)V+=(A[i]-µ)**2
    return V/L
}
,STDEV_P=A=>Math.sqrt(VAR_P(A))
,COVAR_P=(A1,A2)=>{
    const L=A1.length;let i=L,µ1=A1[--i],µ2=A2[i]
    while(--i>=0){µ1+=A1[i];µ2+=A2[i]}µ1/=L;µ2/=L
    i=L;let C=(A1[--i]-µ1)*(A2[i]-µ2)
    while(--i>=0)C+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((C/L).toFixed(3))
}
//Sample Variance S², Standard Deviation S & Covariance
,VAR_S=A=>{
    const L=A.length
    let i=L,µ=A[--i];while(--i>=0)µ+=A[i];µ/=L
    i=L;let V=(A[--i]-µ)**2;while(--i>=0)V+=(A[i]-µ)**2
    return V/(L-1)
}
,STDEV_S=A=>Math.sqrt(VAR_S(A))
,COVAR_S=(A1,A2)=>{
    const L=A1.length;let i=L,µ1=A1[--i],µ2=A2[i]
    while(--i>=0){µ1+=A1[i];µ2+=A2[i]}µ1/=L;µ2/=L
    i=L;let C=(A1[--i]-µ1)*(A2[i]-µ2)
    while(--i>=0)C+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((C/(L-1)).toFixed(3))
}

/*--------------Graphing 📈--------------*/
,TrapeziumRule=(X_Interval,Y_Values_Array)=>Number(((2*Sum(Y_Values_Array)-Y_Values_Array[0]-Y_Values_Array.at(-1))*X_Interval/2).toFixed(3))

,LinearExtrapolation=(x,x0,x1,y0,y1)=>Number((y0+(x-x0)*(y1-y0)/(x1-x0)).toFixed(3))
// Outputs y, also works for extrapolation

,BestFitLine2Points=(x1,x2,y1,y2)=>{
    const m=Number(((y2-y1)/(x2-x1)).toFixed(3))
    return [m,Number((y2-m*x2).toFixed(3))]
    // Outputs [m,c] (m = gradient & c = y-intercep of y=mx+c)
}
,BestFitLineLSM=(X,Y)=>{//Least Square Method: inputs are arrays of x & y coordinates which must be of equal lengths
    let L=X.length,µ_Y=µ_X=N=D=0
    if(L!=Y.length)return'Please enter arrays of x and y coordinates of equal lengths'

    for(let i=L;--i>=0;µ_Y+=Y[i])µ_X+=X[i]
    µ_Y/=L;µ_X/=L

    while(--L>=0){const x=X[L]-µ_X;N+=x*(Y[L]-µ_Y);D+=x**2}
    const m=Number((N/D).toFixed(3))
    return [m,Number((µ_Y-m*µ_X).toFixed(3))]
    // Outputs [m,c] (m = gradient & c = y-intercep of y=mx+c)
}

/*--------------Iterative Greatest Common Divisor & Least Common Multiple--------------*/
,gcd=(a,b)=>{a<b&&([a,b]=[b,a]);while(b!=0)[a,b]=[b,a%b];return a}
,GCD=A=>{let L=A.length,R=A[--L];while(--L>=0)R=gcd(R,A[L]);return R}
,lcm=(a,b)=>a*b/gcd(a,b)
,LCM=A=>{let L=A.length,R=A[--L];while(--L>=0)R=lcm(R,A[L]);return R}

/*--------------Iterative Factorial ! Fractions & Negatives--------------*/
,Gamma=n=>{//The use of this 'Gamma' is to ensure Factorial can work on fractions, not just integers
    // g represents the precision desired, p is the values of p[i] to plug into Lanczos' formula
    //some magic constants
    const g=9,p=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905, -0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7]
    if(n<0.5){
        return Math.PI/Math.sin(n*Math.PI)/Gamma(1-n)
    }else{
        --n;let x=p[0]
        for(let i=0;++i<g;)x+=p[i]/(n+i)
        const t=n+g-1.5
        return Math.sqrt(2*Math.PI)*Math.pow(t,(n+0.5))*Math.exp(-t)*x
    }
}
,Factorial=n=>{
    if(n==0)return 1
    else if(Number.isInteger(n)){
        if(n>0){
            let r=n;while(--n>1)r*=n;return r
        }else//negative integers
            return Infinity
    }else//fraction
        return Gamma(n+1)
}
/*--------------2 Vectors--------------*/
,AngleBetween2Vectors=(A,B)=>{
    let L=A.length,α=A[0],β=B[0],a=α**2,b=β**2,ab=α*β
    if(L!=B.length)return'Please enter A and B arrays of equal lengths'
    while(--L>0){α=A[L];β=B[L];a+=α**2;b+=β**2;ab+=α*β}
    const θ=Math.acos(ab/Math.sqrt(a)/Math.sqrt(b))
    if(θ==0||isNaN(θ))
        console.log('Parallel')
    else if(θ==Math.PI/2)
        console.log('Perpendicular')
    else
        console.log('Not Parallel Or Perpendicular')
    return θ
    // π/2 Radians ≥ Angle ≥ 0
}
,VectorsCrossProduct3D=(A,B)=>[A[1]*B[2]-A[2]*B[1],A[2]*B[0]-A[0]*B[2],A[0]*B[1]-A[1]*B[0]]
// Inputs [A0,A1,A2] & [B0,B1,B2] and returns 3X3 determinant of [[i,j,k],[A0,A1,A2],[B0,B1,B2]] in [i,j,k] format

/*--------------2D Matrices--------------*/
,AddToEachElement=(M,n)=>M.map(r=>r.map(e=>e+=n))
,SubtractFromEachElement=(M,n)=>M.map(r=>r.map(e=>e-=n))
,MultiplyToEachElement=(M,n)=>M.map(r=>r.map(e=>e*=n))
,DivideFromEachElement=(M,n)=>M.map(r=>r.map(e=>e/=n))
,PowerToEachElement=(M,n)=>M.map(r=>r.map(e=>e**=n))

// I can add more of these 👆 but I think you get the general way of making these 😉

,AddMatrices=(M1,M2)=>M1.map((R,r)=>R.map((e,c)=>e+=M2[r][c])) // M1 + M2
,SubtractMatrices=(M1,M2)=>M1.map((R,r)=>R.map((e,c)=>e-=M2[r][c])) // M1 - M2
,DotProductMatrices=(M1,M2)=>{// M1 ⋅ M2 | Matrix Multiplication
    const M2Rows=M2.length,M1Cols=M1[0].length
    if(M1Cols==undefined)//For 1D Arrays SumProduct, correct user input format
        M1=[M1]
    else if(M1[0].length!=M2Rows)
        return'#Cols of 1st Matrix must equal to #rows of 2nd matrix!'

    const Rows=M1.length,Cols=M2[0].length
    ,M=Array.from({length:Rows},()=>Array.from({length:Cols},()=>0))
    // When 2 matrices are multiplied, the resulting matrix has the same # rows as the 1st matrix and the same # cols as 2nd matrix

    // Row x Col → Row , Remember: fiRst seCond thiRd | R.C.R.
    for(let R=-1;++R<Rows;){
        for(let C=-1;++C<Cols;){
            for(let i=-1;++i<M2Rows;){
                M[R][C]+=M1[R][i]*M2[i][C]
            }
        }
    }
    return M
}

,RowElements=(M,r)=>M[r]
,ColElements=(M,c)=>M.map(R=>R[c])
,DiagonalElements=M=>M.map((R,C)=>R[C])

,Determinant2X2=M=>M[0][0]*M[1][1]-M[0][1]*M[1][0] // Only works on Matrix [[a,b],[c,d]] such that det(M)=|M|=a*d-b*c
,Determinant3X3=M=>M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])

,Excluding=i=>xs=>[...xs.slice(0,i),...xs.slice(i+1)]
,RecursiveDeterminant=([xs,...xss])=>xs.length==1?xs[0]:Sum(xs.map((x,i)=>(-1)**i*x*RecursiveDeterminant(xss.map(Excluding(i)))))

// Note: If all the elements in any one row or col in a sqaure matrix are 0, the determinant is 0
,DeterminantOfTriangleMatrix=M=>Number((Product(DiagonalElements(M))).toFixed(3))

,Transpose=M=>M[0].length==1?M.flat():!isNaN(M[0])?M.map(x=>[x]):M[0].map((x,i)=>M.map(x=>x[i]))
// Don't confuse Transpose with Inverse

,IdentityMatrix=i=>Array.from({length:i},(x,r)=>Array.from({length:i},(y,c)=>(r==c?1:0)))

,MinorSubMatrix=(M,Row,Col)=>{// Deletes a row & a col
    if(!Number.isInteger(Row)||!Number.isInteger(Col))return'Please enter Row & Col as integers'
    if(M.length==2)return M[Math.abs(Row-1)][Math.abs(Col-1)]
    M.splice(Row,1);M.forEach(row=>row.splice(Col,1));return M
}

/*--------------Randomness--------------*/
,Shuffle=A=>A.sort(()=>Math.random()-0.5) // Change the order of elements randomly like a deck of cards 🃏🎴

/*--------------Imaginary & Complex Numbers--------------*/
,ComplexSQRT=(a,b)=>{// Square Root √(a+bi), for real a & b
    if([0,undefined].includes(b)){
        return Math.sqrt(a)
    }else{
        const HYPOT_ab=Math.hypot(a,b)
        return [Math.sqrt((HYPOT_ab+a)/2),Math.sign(b)*Math.sqrt((HYPOT_ab-a)/2)] //±[](0) is real, ±[](1) is imaginary
    }
}
,ComplexCBRT=(a,b)=>{// Cube Root ₃√(a+bi), for real a & b
    if([0,undefined].includes(b)){
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
        ,Nature=aNegative?'Maxima':'Minima' // δ²y/δx² = 2a
        if(Discriminant==0){// 1 'repeated' real root which is also the stationary point
            console.log('Root (y=0) & '+Nature+': x = '+x)
        }else{//Either 2 real roots OR 2 complex roots
            const sqrt_discriminantOver_2a=Number((Math.sqrt(Math.abs(Discriminant))/a2).toFixed(3))
            ,y=Number((a*x**2+b*x+c).toFixed(3))
            if(Discriminant>0){//2 real roots
                let x0=Number((x-sqrt_discriminantOver_2a).toFixed(3)),x1=Number((x+sqrt_discriminantOver_2a).toFixed(3))
                aNegative&&([x0,x1]=[x1,x0]) //sort so x₀<x₁
                console.log('• Roots (y=0): x₀ = '+x0+' , x₁ = '+x1+'\n• '+Nature+': x = '+x+' , y = '+y)
            }else{//2 complex roots
                let Roots=x+' ± '+Math.abs(sqrt_discriminantOver_2a)+' i';Math.abs(sqrt_discriminantOver_2a)==1&&(Roots=x+' ± i')
                console.log('• Roots (y=0): x = '+Roots+'\n• '+Nature+': x = '+x+' , y = '+y)
            }
        }
    }
}
,QuadraticDiscriminant=(a,b,c)=>{//y=ax²+bx+c, also known as parabolic
    if(a==0){
        return SolveLinear(b,c)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_2
        const Discriminant=b**2-4*a*c
        if(Discriminant==0){
            console.log("• 1 'repeated' real root which is also the stationary point")
        }else if(Discriminant>0){
            console.log('• 2 real roots')
        }else{//Discriminant<0
            console.log('• 2 complex roots')
        }
        return Discriminant
    }
}
,CubicDiscriminant=(a,b,c,d)=>{//y=ax³+bx²+cx+d
    if(a==0){
        return SolveQuadratic(b,c,d)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_3
        const Discriminant=18*a*b*c*d-4*d*b**3+b**2*c**2-4*a*c**3-27*a**2*d**2
        if(Discriminant==0){
            console.log('• At least 1 stationary point is also a root')
        }else if(Discriminant>0){//∴ Quadratic_Discriminant 4*(b**2-3*a*c)>0
            console.log('• 3 real distinct roots')
        }else{//Discriminant<0
            console.log("• Only 1 real root which isn't a stationary point")
        }
        return Discriminant
    }
}
,QuarticDiscriminant=(a,b,c,d,e)=>{//y=ax⁴+bx³+cx²+dx+e
    if(a==0){
        return CubicDiscriminant(b,c,d,e)
    }else{//https://en.wikipedia.org/wiki/Discriminant#Degree_4
        const Discriminant=256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*c**2*e**2+144*a**2*c*d**2*e-27*a**2*d**4+144*a*b**2*c*e**2-6*a*b**2*d**2*e-80*a*b*c**2*d*e+18*a*b*c*d**3+16*a*c**4*e-4*a*c**3*d**2-27*b**4*e**2+18*b**3*c*d*e-4*b**3*d**3-4*b**2*c**3*e+b**2*c**2*d**2
        if(Discriminant==0){
            console.log('• At least 1 stationary point is also a root')
        }else if(Discriminant>0){
            console.log('• Roots are either all real or all complex')
        }else{//Discriminant<0
            console.log('• 2 real roots and 2 complex roots')
        }
        return Discriminant
    }
}

//Note: It has been proven in 1824 in the Abel–Ruffini theorem that there cannot be a general solution for polynomials of degrees greater than 4

/*--------------Logarithms & Exponentials--------------*/ // log(😅) = 💧log(😄)
,BaseLog=(Base,Num)=>Number((Math.log(Num)/Math.log(Base)).toFixed(3))

/*--------------Convert Bases--------------*/
/*
    There are few ways to convert between bases
    The fastest way to convert, which also never result in non-terminating repeating decimals involves replacing each character with others,
    and therefore doesn't involve any actual math calculations.

    For 2 bases to quality to this efficient conversion, their BaseLog needs to be an integer.
    Their BaseLog integer specifies the number of characters which are replaced in the smaller base from a single character in the larger base.
    2 Bases can also use this method if thier BaseLog isn't an integer, but they share another base to which they seperately have integer BaseLog,
    e.g. Base8 Octal and Base16 Hexadecimal are efficiently conveted between via Base2 Binary.
*/

,Base4ToBase2=InputNumber=>{
    InputNumber=String(InputNumber);let i=-1,Point=false,Result=''
    if(InputNumber[0]=='-'){InputNumber=InputNumber.slice(1);Result='-'}
    const L_1=InputNumber.length-1
    do{
        const Digit=InputNumber[++i]
        switch(Digit){// BaseLog(2,4)=2, hence each digit in Base4 is replaced with exactly 2 digits in Base2
            case '0':Result+='00';break
            case '1':Result+='01';break
            case '2':Result+='10';break
            case '3':Result+='11';break
            case '-':return"Can only have a minus sign '-' at the start"
            case '.':if(Point){return"Can't have multiple decimal points '.'"}else{Result+='.';Point=true;break}
            default: return"Character '"+Digit+"' is invalid"
        }
    }while(i<L_1)
    return Result
}
,Base8ToBase2=InputNumber=>{
    InputNumber=String(InputNumber);let i=-1,Point=false,Result=''
    if(InputNumber[0]=='-'){InputNumber=InputNumber.slice(1);Result='-'}
    const L_1=InputNumber.length-1
    do{
        const Digit=InputNumber[++i]
        switch(Digit){// BaseLog(2,8)=3, hence each digit in Base8 is replaced with exactly 3 digits in Base2
            case '0':Result+='000';break
            case '1':Result+='001';break
            case '2':Result+='010';break
            case '3':Result+='011';break
            case '4':Result+='100';break
            case '5':Result+='101';break
            case '6':Result+='110';break
            case '7':Result+='111';break
            case '-':return"Can only have a minus sign '-' at the start"
            case '.':if(Point){return"Can't have multiple decimal points '.'"}else{Result+='.';Point=true;break}
            default: return"Character '"+Digit+"' is invalid"
        }
    }while(i<L_1)
    return Result
}
// This allows converting both integers and decimals from any base to any base in the range of 2-36
,ConvertBases=(InputNumber,InputBase,OutputBase)=>{
    if(!Number.isInteger(InputBase)||!Number.isInteger(OutputBase)||InputBase>36||InputBase<2||OutputBase>36||OutputBase<2){
        return'InputBase & OutputBase must be whole numbers between 2-36'
    }else if(InputBase==OutputBase){
        console.log('InputBase = OutputBase');return InputNumber
    }else{
        const AllowedCharcters=['.','-','0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        ,Base10=InputBase==10?InputNumber:parseInt(InputNumber,InputBase)
        ,StringInput=String(InputNumber)

        InputBase==36||AllowedCharcters.splice(InputBase-36)
        let L=StringInput.length
        do{
            const Char=StringInput[--L]
            if(!AllowedCharcters.includes(Char))return"Character '"+Char+"' is invalid in base "+InputBase
        }while(L>0)

        const Result=Base10.toString(OutputBase),ResultNumber=Number(Result)
        return isNaN(ResultNumber)?Result:ResultNumber
    }
}

/*--------------Convert Units--------------*/
,ConvertMassUnits=(InputNumber,UnitIndex)=>{ // UnitIndex would be one of the keys from the Units object such as 'kg' or 'lb'
    UnitIndex==undefined&&(UnitIndex='g') // Gram is default if no unit is specified
    const Units={
        // Mass UnitIndices:

        // British Imperial:
        'gr': 480/31.1034768 // Grain
        ,'dwt': 20/31.1034768 // Pennyweight
        ,'ozt': 1/31.1034768 // Troy Ounce

        // Metric:
        ,'ng': 10**9 // Nanogram
        ,'mcg': 10**6 // Microgram , A.K.A. 'μg', but 'μ' isn't used here because non-UTF-8 characters sometimes render badly in some consoles
        ,'mg': 1000 // Milligram
        ,'ct': 5 // Carat
        ,'g': 1 // Gram - reference
        ,'kg': 0.001 // Kilogram
        ,'t': 10**-6 // Tonne

        // American Imperial:
        ,'oz': 16/453.59237 // Ounce
        ,'lb': 1/453.59237 // Pound
        ,'st': 1/14/453.59237 // Stone
        ,'tn': 7/14000/453.59237 // US Ton
        ,'LT': 7/14000/453.59237/1.12 // Imperial Long Ton
    }
    ,Ratio=InputNumber/Units[UnitIndex];Units[UnitIndex]=InputNumber
    Object.keys(Units).forEach(Key=>Key==UnitIndex||(Units[Key]=Number((Units[Key]*Ratio).toPrecision(4))))
    return Units
}

/*--------------Standard Maclaurin Series--------------*/
// Note: The functions below all have an equivilant Math API function built into JS, but I consider these as syntactic sugar and redundant, so here's what they would look like if we only had the JS syntax we acutally NEEDED for basic functionality

// x = number , i = # iterations (more takes more time, but increases accuracy)

,Exponential=(x,i)=>{// Math.exp(x) | Math.E**x | e^x
    if(i==undefined)
        i=30 // default # iterations
    else if(!Number.isInteger(i)||i<1)
        return'i must be a positive integer!'

    let j=1,N=x,r=x+1
    while(++j<=i){N*=x/j;r+=N}
    return r
}

/*--------------Strings--------------*/
,FindIndices=(SearchedElement,A)=>{
    const Indices=[],L=A.length;let i=-1;do{A[i]==SearchedElement&&Indices.push(i)}while(++i<L);return Indices
} // e.g. FindIndices('h',[true,'h',6,'try',6,false,'h',0,-5.4,'h']) => [1, 6, 9]

,RemoveLeadingCharacter=(Str,Char)=>{while(Str[0]==Char)Str=Str.slice(1);return Str}
,RemoveFollowUpCharacter=(Str,Char)=>{while(Str.at(-1)==Char)Str=Str.slice(0,-1);return Str}

```
