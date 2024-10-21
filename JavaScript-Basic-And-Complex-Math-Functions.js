'use strict'
/*--------------Statistics--------------*/
// Key: A = Array, i.e. [#,#,..]

const Minimum=A=>{let L=A.length,M=A[--L];while(--L>=0){const N=A[L];N<M&&(M=N)};return M}
,Maximum=A=>{let L=A.length,M=A[--L];while(--L>=0){const N=A[L];N>M&&(M=N)};return M}
,MinMax3N2=A=>{//This is a slightly faster algorithm to get both the Minimum and Maximum of an array at the same time
    let L=A.length,Min,Max
    if(L&1)//if odd
        Min=Max=A[--L]
    else{//if even
        const a=A[--L],b=A[--L];a>b?(Max=a,Min=b):(Max=b,Min=a)
    }
    for(let Big,Small;L>0;){
        const a=A[--L],b=A[--L]
        a>b?(Big=a,Small=b):(Big=b,Small=a)
        Big>Max&&(Max=Big);Small<Min&&(Min=Small)
    }
    return[Min,Max]
}
,Range=A=>{const a=MinMax3N2(A);return a[1]-a[0]} // Max-Min

,Sum=A=>{let L=A.length,S=0;do{S+=A[--L]}while(L>0);return S} /* Σ
  e.g. Sum([3,4,4,5,12]) ➜ 28 */

,MeanAverage=A=>{const L=A.length;let l=L,Σ=A[--l];while(--l>=0)Σ+=A[l];return Σ/L} /* µ
   e.g. MeanAverage([3,4,5,12]) ➜ 6 , MeanAverage([3,4,4,5,12]) ➜ 5.6 */

,MedianAverage=A=>{A.sort((a,b)=>a-b);const L2=A.length/2;return L2%1?A[L2-0.5]:(A[L2]+A[L2-1])/2}
// e.g. MedianAverage([3,4,5,12]) ➜ 4.5 , MedianAverage([3,4,4,5,12]) ➜ 4

,Unique=A=>[...new Set(A)] /* e.g. Unique([3,4,'4',5,4,12]) ➜ [3,4,'4',5,12]
    Alternatively: Unique=A=>Array.from(new Set(A)) OR Unique=A=>A.filter((e,i,s)=>s.indexOf(e)==i) */
,CountOccurances=A=>{const C={};for(let L=A.length;--L>=0;){const E=A[L];C[E]?C[E]+=1:C[E]=1}return C}
// e.g. CountOccurances([3,4,4,5,12]) ➜ {'3':1,'4':2,'5':1,'12':1}

,ModeAverage=A=>{//Array items with higest occurances, not just numbers but also strings
    Array.isArray(A)&&(A=CountOccurances(A))//Input can be array [] or object {}
    const HighestOccurance=Maximum(Object.values(A))
    ,keys=Object.keys(A)
    let Mode=filteredKeys=keys.filter(key=>A[key]==HighestOccurance)

    Mode.length<2&&(Mode=Mode[0])/* if there's no tie for the highest occurance, 
    return only that value, otherwise return an array of the tied values */
    // May Do: convert each tied mode value to number if numeric & console.log 'tie' to inform

    console.log('Mode: ',Mode,', Occurance: '+HighestOccurance)
    return[Mode,HighestOccurance]
}// e.g. ModeAverage([3,4,4,5,12]) ➜ 'Mode:  4 , Occurance: 2' , [ '4', 2 ]

,Product=A=>{// ∏
    let L=A.length,P=A[--L];if(P==0)return 0
    while(--L>=0){const a=A[L];if(a==0)return 0;P*=a}
    return P
}// e.g. Product([1,3,4,7]) ➜ 84 , 1*3*4*7

// Efficiently returns 0 as soon as any is detected. This can also be modified to detect Infinity, -Infinity and NaN

,SumProduct=(A1,A2)=>{
    let L=A1.length,R=A1[--L]*A2[L]
    while(--L>=0)R+=A1[L]*A2[L]
    return R
}/* e.g. SumProduct([1,3],[4,7]) ➜ 25 , 1*4+3*7
    Weight Sum Model (WSM) in Multi-Criteria Analysis (MCA)
*/
,SumPower=(A1,A2)=>{
    let L=A1.length,R=A1[--L]**A2[L]
    while(--L>=0)R+=A1[L]**A2[L]
    return R
}// e.g. SumPower([1,3],[4,7]) ➜ 2188 , 1**4+3**7

,ProductPower=(A1,A2)=>{
    let L=A1.length,R=A1[--L]**A2[L];if(R==0)return 0
    while(--L>=0){const a=A1[L];if(a==0)return 0;R*=a**A2[L]}
    return R
}/* e.g. ProductPower([1,3],[4,7]) ➜ 2187 , 1**4*3**7
    Weight Product Model (WPM) in Multi-Criteria Analysis (MCA)
*/

,AscendingSort=A=>A.sort((a,b)=>a<b?-1:1) /* Both Numbers and Strings Alphabetically Ascendingly
    e.g. AscendingSort([-7.6,0,0.5,1,4.3,true,false,NaN,'str',' ',[],{},Infinity,-Infinity]) =>
        [-Infinity,-7.6,0,false,[],' ',0.5,1,true,4.3,NaN,{},'str',Infinity] */
,SortNumbers=A=>A.sort((a,b)=>a-b) /* Only numbers, but slightly more efficent. From smallest to biggest
    e.g. SortNumbers([0,-1.4,7.6]) ➜ [-1.4,0,7.6] */

,Rank=A=>{const S=A.slice().sort((a,b)=>a>b?-1:1);return A.map(v=>S.indexOf(v)+1)}
// rank values with biggest being #1 , e.g. Rank([0,-1.4,7.6]) ➜ [2,3,1] , Rank(['a','c','b']) ➜ [3,1,2]

,FractionRank=A=>{// Adds 0.5*(occurrences-1) per each tied rank
    const R=Rank(A),L=R.length,C=new Map()
    for(let i=-1;++i<L;){const r=R[i];C.has(r)||C.set(r,0);C.set(r,C.get(r)+1)}
    for(let i=-1;++i<L;)R[i]+=(C.get(R[i])-1)/2
    return R
}//e.g. FractionRank([5,62,5,0,4,-3.4,5,62]) ➜ [4, 1.5, 4, 7, 6, 8, 4, 1.5]

,SpearmanCorrelation=(A1,A2)=>{// Coefficient 'ρ': 1 ≥ ρ ≥ -1
    const n=A1.length,Rank1=FractionRank(A1),Rank2=FractionRank(A2);let Σd2=0
    for(let L=n;--L>=0;)Σd2+=(Rank1[L]-Rank2[L])**2
    return Number((1-(6*Σd2/(n*(n**2-1)))).toFixed(3)) // ρ = 1 - 6Σd²/n(n²-1)
}// e.g. SpearmanCorrelation([1,6,4,3,4],[9,5,9,7,2]) ➜ -0.475

,PearsonCorrelation=(A1,A2)=>{// Coefficient 'r': 1 ≥ r ≥ -1
    const n=A1.length
    if(n<2)return -Infinity // to prevent error in the console when used in some apps that may not accept NaN
    let i=n,S1=A1[--i],S2=A2[i],SP11=S1**2,SP22=S2**2,SP12=S1*S2 // 'S' is Sum and 'SP' is SumProduct
    while(--i>=0){
        const a1=A1[i],a2=A2[i]
        S1+=a1;S2+=a2;SP11+=a1**2;SP22+=a2**2;SP12+=a1*a2
    }
    return(n*SP12-S1*S2)/Math.sqrt((n*SP11-S1**2)*(n*SP22-S2**2))
    // r = (n∙Σ(x∙y)-Σ(x)∙Σ(y))/√((n∙Σ(x²)-(Σx)²)∙(n∙Σ(y²)-(Σy)²))
}// e.g. PearsonCorrelation([1,6,3,4],[9,5,7,2]) ➜ -0.520

//Population Covariance, Variance σ² & Standard Deviation σ
,COVAR_P=(A1,A2)=>{
    const L=A1.length;let i=L,µ1=A1[--i],µ2=A2[i]
    while(--i>=0){µ1+=A1[i];µ2+=A2[i]}µ1/=L;µ2/=L
    i=L;let C=(A1[--i]-µ1)*(A2[i]-µ2)
    while(--i>=0)C+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((C/L).toFixed(3))
}
,VAR_P=A=>{
    const L=A.length
    let i=L,µ=A[--i];while(--i>=0)µ+=A[i];µ/=L
    i=L;let V=(A[--i]-µ)**2;while(--i>=0)V+=(A[i]-µ)**2
    return Number((V/L).toFixed(3))
}/* e.g. VAR_P([85,92,64,99,56]) ➜ 271.76
    source: https://study.com/skill/learn/calculating-population-standard-deviation-explanation.html */
,STDEV_P=A=>Number(Math.sqrt(VAR_P(A)).toFixed(3))// e.g. STDEV_P([85,92,64,99,56]) ➜ 16.485

//Sample Covariance, Variance S² & Standard Deviation S
,COVAR_S=(A1,A2)=>{
    const L=A1.length;let i=L,µ1=A1[--i],µ2=A2[i]
    while(--i>=0){µ1+=A1[i];µ2+=A2[i]}µ1/=L;µ2/=L
    i=L;let C=(A1[--i]-µ1)*(A2[i]-µ2)
    while(--i>=0)C+=(A1[i]-µ1)*(A2[i]-µ2)
    return Number((C/(L-1)).toFixed(3))
}
,VAR_S=A=>{
    const L=A.length
    let i=L,µ=A[--i];while(--i>=0)µ+=A[i];µ/=L
    i=L;let V=(A[--i]-µ)**2;while(--i>=0)V+=(A[i]-µ)**2
    return Number((V/(L-1)).toFixed(3))
}// e.g. VAR_S([9,2,5,4,12,7]) ➜ 13.1 , source: https://www.ztable.net/sample-population-standard-deviation
,STDEV_S=A=>Number(Math.sqrt(VAR_S(A)).toFixed(3))// e.g. STDEV_S([9,2,5,4,12,7]) ➜ 3.619

/* Use Sample S equations of you only analyse a subset of a larger population, and Population σ equations if you analyse the entire population.
Variance returns the same result as entiring the same array as the 2 inputs for Covariance */

,FindIndices=(SearchedElement,A)=>{
    const Indices=[],L=A.length;let i=-1;do{A[i]==SearchedElement&&Indices.push(i)}while(++i<L);return Indices
} // e.g. FindIndices('h',[true,'h',6,'try',6,false,'h',0,-5.4,'h']) ➜ [1, 6, 9]

/*--------------Multi-Criteria Analysis (MCA)--------------*/

,NormaliseWeights=Weights=>{
    // This function inputs an array of numbers and divides each by their total, so the returned array will have a sum of 1.00

    let L=l=Weights.length,S=0;do{S+=Weights[--L]}while(L>0)
    while(--l>=0)Weights[l]/=S
    return Weights
}// e.g. NormaliseWeights([3,3,1,1,1,1,2]) ➜ [0.25,0.25,0.083,0.083,0.083,0.083,0.166]
,NormaliseScores=(Scores,Best,Worst)=>{
    // For a given single criterion 'Scores' is an array where each value corresponds to an alternative

    /*----------------------Checks----------------------*/ 
    if(isNaN(Best))return'Best Score Must Be A Number!'
    if(isNaN(Worst))return'Worst Score Must Be A Number!'
    if(Best==Worst)return'Best And Worst Scores Must Be Different Numbers!'

    /*----------------------Calculations----------------------*/ 
    for(let L=Scores.length;--L>=0;){
        const Score=(Scores[L]-Worst)/(Best-Worst) // Normalise the score between 1 and 0

        Scores[L]=Score<=0?0:Score>1?1:Score-0/* 
            if a score is out of bounds, make it equal to its nearest bound, i.e. 1 OR 0 
            '<=0' is used to avoid a score of '-0'
            'Score-0' is to return NaN instead of scripting isNaN(Score) to minify 
        */
    }
    return Scores
}// e.g. NormaliseScores([3,7,8,2,5,6,'String'],7,3) ➜ [0,1,1,0,0.5,0.75,NaN]

,PROMETHEE=(W,S,II)=>{
    // II: if(II==true){PROMETHEE II}else{PROMETHEE I}; W = Array of weights, S = Matrix of normalised scores between 0 and 1, with 0 being worst and 1 being best. # Rows = # Alternatives , # Cols = # Criteria Weights.

    /*----------------------Checks----------------------*/ 
    if(!W.every(s=>s>=0))return'Every weight needs to be either a positive number or 0'
    if(!S.flat(Infinity).every(s=>s<=1&&s>=0))
        return'Every score needs to be pre-normalised between 0 and 1. Please consider using the NormaliseScores function.'
    const C=W.length // C = # Criteria
    if(S[0].length!=C)return'#Criteria must equal the #Weights'

    /*----------------------Preparations----------------------*/ 
    const A=S.length // A = # Alternatives
    ,F=Array.from({length:3},()=>Array(A).fill(0)) // F has 3 cols for flows: Leaving, Entering And Net Total
    II=II==true?A-1:1

    /*----------------------Calculations----------------------*/
    for(let Am=A;--Am>=0;)// Am = Minuend Alternative
        for(let As=A;--As>=0;)// As = Subtrahend Alternative
            if(Am!=As)
                for(let c=C;--c>=0;){// c = Column
                    const Flow=Math.max(0,W[c]*(S[Am][c]-S[As][c]))/II
                    F[0][Am]+=Flow;F[1][As]+=Flow;F[2][Am]+=Flow;F[2][As]-=Flow
                }
    F[3]=Rank(F[2]) // last col is ranks with #1 being the best
    return F // returns a matrix of 4 columns with # Rows = # Alternatives
}
// e.g. PROMETHEE([0.38,0.09,0.05,0.22,0.17,0.08,0.02],[[0,0,0,0,0,1,0],[0.5,0.5,1/3,0,2/3,0.14286,0],[0.5,1,1,0.5,2/3,0.71429,1/3],[1,1,1,1,1,0,1],[0.5,0.5,1,0.5,2/3,1,1],[0,0,1/3,0,2/3,0.42857,0]],true) - https://www.econstor.eu/bitstream/10419/237008/1/1752580664.pdf

/*--------------Generate 1D Arrays--------------*/
/*
To generate an array of number 
The functions below are inspired by Casio calculators:
*/

,GenerateArraySSE=(Start,Step,End)=>{
    const A=[]
    do{A.push(Start);Start+=Step}while(Start<=End)
    return A
}// e.g. GenerateArraySSE(1,1.5,9) ➜ [1,2.5 4,5.5,7,8.5]

,GenerateArraySLS=(Start,Length,Step)=>{
    const A=[Start]
    for(let i=0;++i<Length;A[i]=Start)Start+=Step
    return A
}// e.g. GenerateArraySLS(0,10,1) ➜ [0,1,2,3,4,5,6,7,8,9]

,GenerateArrayLSE=(Length,Step,End)=>{
    const A=[]
    while(--Length>=0){A[Length]=End;End-=Step}
    return A
}// e.g. GenerateArrayLSE(8,2,30) ➜ [16,18,20,22,24,26,28,30]

/*--------------1D Arrays & 2D Matrices--------------*/

,ArrayFuncForEachElement=(F,A)=>A.map(x=>F(x)) // apply function F for each element x in array A
,MatrixFuncForEachElement=(F,M)=>M.map(r=>r.map(x=>F(x))) // apply function F for each element x in matrix M

,ArrayAddToEachElement=(A,n)=>A.map(e=>n+e) //e.g. ArrayAddToEachElement([-4,0,5.5],3) ➜ [-1,3,8.5]
,MatrixAddToEachElement=(M,n)=>M.map(r=>r.map(e=>e+=n)) //e.g. ArrayAddToEachElement([[-4,0],[5.5,-2.1]],3) ➜ [-1,3,8.5]

,ArraySubtractFromEachElement=(A,n)=>A.map(e=>e-n) //e.g. ArraySubtractFromEachElement([-4,0,5.5,-2.1],3) ➜ [-7,-3,2.5,-5.1]
,MatrixSubtractFromEachElement=(M,n)=>M.map(r=>r.map(e=>e-n)) //e.g. MatrixSubtractFromEachElement([[-4,0],[5.5,-2.1]],3) ➜ [[-7,-3],[2.5,-5.1]]

,ArrayMultiplyToEachElement=(A,n)=>A.map(e=>e*n) //e.g. ArrayMultiplyToEachElement([-4,0,5.5,-2.1],3) ➜ [-12,0,16.5,-6.3]
,MatrixMultiplyToEachElement=(M,n)=>M.map(r=>r.map(e=>e*n)) //e.g. MatrixMultiplyToEachElement([[-4,0],[5.5,-2.1]],3) ➜ [[-12,0],[16.5,-6.3]]

// I can add more of these 👆 but I think you get the general way of making these 😉
// Remember: you can turn any matrix into an array using A.flat() or A.flat(Infinity)

/*--------------2D Matrices--------------*/

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

,Transpose=M=>M[0].length==1?M.flat():!isNaN(M[0])?M.map(x=>[x]):M[0].map((x,i)=>M.map(x=>x[i]))
// Don't confuse Transpose with Inverse

,IdentityMatrix=i=>Array.from({length:i},(x,r)=>Array.from({length:i},(y,c)=>(r==c?1:0)))

,MinorSubMatrix=(M,Row,Col)=>{// Deletes a row & a col
    if(!Number.isInteger(Row)||!Number.isInteger(Col))return'Please enter Row & Col as integers'
    if(M.length==2)return M[Math.abs(Row-1)][Math.abs(Col-1)]
    M.splice(Row,1);M.forEach(row=>row.splice(Col,1));return M
}

/*
Only Square (#Rows=#Cols) Matrices with non-zero determinant have an inverse matrix. This is because M⁻¹ = 1/det(M) ⋅ C, where C is the cofactor matrix, and only square matrices have determinants, if the determinant is zero, it would cause a division by zero as the determinant is a donominator.

An inverse of a matrix is defined as: M ⋅ M⁻¹ =  M⁻¹ ⋅ M = I , where I is an identity matrix of the same length.
Consequently, and as an aside, the inverse of the inverse matrix is equal to the original matrix: (M⁻¹)⁻¹ = M
*/

/*--------------Graphing 📈--------------*/
,TrapeziumRule=(x,Y)=>Number(((2*Sum(Y)-Y[0]-Y.at(-1))*x/2).toFixed(3)) // x = x-interval , Y = y-values Array

,LinearPolation=(x,x0,x1,y0,y1)=>{
    x>x0&&x>x1||x<x0&&x<x1?console.log('Extrapolation'):console.log('Interpolation')
    return Number((y0+(x-x0)*(y1-y0)/(x1-x0)).toFixed(3))// Outputs y
}// e.g. LinearPolation(5,1,8,3,8) ➜ Interpolation,5.857 , LinearPolation(-5,1,8,-3,8) ➜ Extrapolation,-12.429

,BestFitLine2Points=(x1,x2,y1,y2)=>{// Outputs [m,c] (m = gradient slope & c = y-intercep, of y = m⋅x + c)
    const m=Number(((y2-y1)/(x2-x1)).toFixed(3))
    return[m,Number((y1-m*x1).toFixed(3))] // alternatively (y2-m*x2)
}/* e.g. BestFitLine2Points(2,6,3,5) ➜ [0.5,2] , i.e. y = x/2 + 2
Other useful linear:
    m = (y - c) / x , for a specific single (x,y) coordinate with known c
    c = y - m⋅x , for a specific single (x,y) coordinate with known gradient
    x-intercep = -c/m , once both m and c are known
    e.g. -2/0.5 ➜ -4 
*/
,BestFitLineLSM=(X,Y)=>{//Least Square Method: inputs are arrays of x & y coordinates which must be of equal lengths
    let L=X.length,µ_Y=µ_X=N=D=0
    if(L!=Y.length)return'Please enter arrays of x and y coordinates of equal lengths'

    for(let i=L;--i>=0;µ_Y+=Y[i])µ_X+=X[i]
    µ_Y/=L;µ_X/=L

    while(--L>=0){const x=X[L]-µ_X;N+=x*(Y[L]-µ_Y);D+=x**2}
    const m=Number((N/D).toFixed(3))
    return[m,Number((µ_Y-m*µ_X).toFixed(3))]
}// Outputs [m,c] (m = gradient slope & c = y-intercep, of y = m⋅x + c)

/*--------------Iterative Greatest Common Divisor & Least Common Multiple--------------*/
/*
    a & b are single numbers, use gcd and lcm when only dealing with 2 numbers
    A is an array of numbers, use GCD and LCM when dealing with 3 or more numbers
*/
,gcd=(a,b)=>{let c;if(a<b){c=a;a=b;b=c};while(b!=0){c=a;a=b;b=c%b};return a} // e.g. gcd(24,32) ➜ 8
,GCD=A=>{let L=A.length,R=A[--L];while(--L>=0)R=gcd(R,A[L]);return R} // e.g. GCD([16,40,64]) ➜ 8
,lcm=(a,b)=>a*b/gcd(a,b) // e.g. lcm(6,9) ➜ 18
,LCM=A=>{let L=A.length,R=A[--L];while(--L>=0)R=lcm(R,A[L]);return R} // e.g. LCM([6,9,12]) ➜ 36

/*--------------Iterative Factorial ! Fractions & Negatives--------------*/
,Gamma=n=>{//The use of this 'Gamma' is to ensure Factorial can work on fractions, not just integers
    // g represents the precision desired, p is the values of p[i] to plug into Lanczos' formula
    //some magic constants
    const g=9,p=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905, -0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7]
    if(n<0.5)return Math.PI/Math.sin(n*Math.PI)/Gamma(1-n)
    else{
        --n;let x=p[0]
        for(let i=0;++i<g;)x+=p[i]/(n+i)
        const t=n+g-1.5
        return Math.sqrt(2*Math.PI)*Math.pow(t,(n+0.5))*Math.exp(-t)*x
    }
}
,Factorial=n=>{
    if(n==0)return 1
    else if(Number.isInteger(n)){
        if(n>0){//positive integers
            let r=n;while(--n>1)r*=n;return r
        }else//negative integers
            return Infinity
    }else//fraction
        return Gamma(n+1)

/*--------------Iterative Fibonacci & Golden Ratio φ--------------*/
,Fib=I=>{//Fibonacci , I = Integer > -1
    if(I==0)return 0;if(I==1)return 1
    if(!Number.isInteger(I)||I<0)return'Please Enter A Positive Integer!'
    let a=0,b=i=1,f;while(++i<=I){f=a+b;a=b;b=f}
    return f
}// e.g. Fib(7) ➜ 13 , Fib(8) ➜ 21 , Fib(9) ➜ 34

,φ=(1+5**0.5)/2 /* Golden Ratio 1.618033988749895..
    Can be derived as a ratio of Fib(I+1)/Fib(I) where I is a large integer
    e.g. Fib(41)/Fib(40) ➜ 1.618033988749895
    I found by trial & error that I = 40 is the smallest integer to yield this accuracy
*/

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

/*--------------Randomness--------------*/
,Shuffle=A=>A.sort(()=>Math.random()-0.5) // Change the order of elements randomly like a deck of cards 🃏🎴

/*--------------Imaginary & Complex Numbers--------------*/
,ComplexSQRT=(a,b)=>{// Square Root √(a+bi), for real a & b
    if(b==undefined) // if no b is entered
        return Math.sqrt(a)
    // else if a non-zero b number is entered
        const HYPOT_ab=Math.hypot(a,b)
        return[Math.sqrt((HYPOT_ab+a)/2),Math.sign(b)*Math.sqrt((HYPOT_ab-a)/2)]
        //±[,][0] is real, ±[,][1] is imaginary
}
,ComplexCBRT=(a,b)=>{// Cube Root ₃√(a+bi), for real a & b
    if(b==undefined) // if no b is entered
        return Math.cbrt(a)
    // else if a non-zero b number is entered
        const θ=Math.atan2(b,a)/3,r=Math.cbrt(Math.hypot(a,b))
        return[r*Math.cos(θ),r*Math.sin(θ)] //[,][0] is real, [,][1] is imaginary
}

/*---Polynomials: General Solutions To Find Roots, Stationary Points & Their Natures---*/
//These find all the real and complex roots for x @ y=0 and stationary points and their natures for real/non-complex coefficients a, b, c, d, e & f

,SolveLinear=(m,c)=>//y=mx+c , m = gradient, c = y-intercep
    m==0?console.log('No Gradient? No Root!'):console.log('Root (y=0): x = '+Number((-c/m).toFixed(3)))

,SolveQuadratic=(a,b,c)=>{//y=ax²+bx+c, also known as parabolic
    if(a==0)//Must use SolveLinear(b,c) @ a=0, otherwise you'd divide by 0
        return SolveLinear(b,c)
    else{//https://en.wikipedia.org/wiki/Discriminant#Degree_2
        const Discriminant=b**2-4*a*c,a2=2*a,aNegative=a<0,x=Number((-b/a2).toFixed(3))// δy/δx = 0 = 2ax+b ∴ x=-b/(2a)
        ,Nature=aNegative?'Maxima':'Minima' // δ²y/δx² = 2a
        if(Discriminant==0)// 1 'repeated' real root which is also the stationary point
            console.log('Root (y=0) & '+Nature+': x = '+x)
        else{//Either 2 real roots OR 2 complex roots
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
    if(a==0)
        return SolveLinear(b,c)
    else{//https://en.wikipedia.org/wiki/Discriminant#Degree_2
        const Discriminant=b**2-4*a*c
        if(Discriminant==0)
            console.log("• 1 'repeated' real root which is also the stationary point")
        else if(Discriminant>0)
            console.log('• 2 real roots')
        else//Discriminant<0
            console.log('• 2 complex roots')
        return Discriminant
    }
}
,CubicDiscriminant=(a,b,c,d)=>{//y=ax³+bx²+cx+d
    if(a==0)
        return SolveQuadratic(b,c,d)
    else{//https://en.wikipedia.org/wiki/Discriminant#Degree_3
        const Discriminant=18*a*b*c*d-4*d*b**3+b**2*c**2-4*a*c**3-27*a**2*d**2
        if(Discriminant==0)
            console.log('• At least 1 stationary point is also a root')
        else if(Discriminant>0)//∴ Quadratic_Discriminant 4*(b**2-3*a*c)>0
            console.log('• 3 real distinct roots')
        else//Discriminant<0
            console.log("• Only 1 real root which isn't a stationary point")
        return Discriminant
    }
}
,QuarticDiscriminant=(a,b,c,d,e)=>{//y=ax⁴+bx³+cx²+dx+e
    if(a==0)
        return CubicDiscriminant(b,c,d,e)
    else{//https://en.wikipedia.org/wiki/Discriminant#Degree_4
        const Discriminant=256*a**3*e**3-192*a**2*b*d*e**2-128*a**2*c**2*e**2+144*a**2*c*d**2*e-27*a**2*d**4+144*a*b**2*c*e**2-6*a*b**2*d**2*e-80*a*b*c**2*d*e+18*a*b*c*d**3+16*a*c**4*e-4*a*c**3*d**2-27*b**4*e**2+18*b**3*c*d*e-4*b**3*d**3-4*b**2*c**3*e+b**2*c**2*d**2
        if(Discriminant==0)
            console.log('• At least 1 stationary point is also a root')
        else if(Discriminant>0)
            console.log('• Roots are either all real or all complex')
        else//Discriminant<0
            console.log('• 2 real roots and 2 complex roots')
        return Discriminant
    }
}

//Note: It has been proven in 1824 in the Abel–Ruffini theorem that there cannot be a general solution for polynomials of degrees greater than 4

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

/*--------------Fraction--------------*/
,DecimalToFraction=d=>{// d = Decimal
    const ε=1e-10,S=Math.sign(d)// S = sign, meaning -1, 0 or 1
    let ABS=Math.abs(d)
    const G=ABS>1;G&&(ABS**=-1)
    if(isNaN(S))return'Please input a number'
    let a=Math.floor(ABS),i=h2=N=0,D=k2=1,y=ABS-a,x,h,k

    while(y>ε){
        if(++i==100)return'100 iterations passed without finding a ratio of integers, so the input number is probably irrational'
        x=1/y;a=Math.floor(x);y=x-a
        h=a*D+h2;h2=D;D=h
        k=a*N+k2;k2=N;N=k
    }
    console.log(i,' iterations');return G!=true?[S*N,D]:[S*D,N]// [ N = Numerator , D = Denominator ]
}/*e.g. DecimalToFraction(76.92) ➜ [1923,25]; DecimalToFraction(-0.0314) ➜ [-157,5000]
    The uniquely nice thing about this function is that the output is equal to the input,so you can easily verify
    the outputs are correct by entering them into a terminal to get the input, e.g. 1923/25 ➜ 76.92 , -157/5000 ➜ -0.0314
*/

/*--------------Logarithms & Exponentials--------------*/ // log(😅) = 💧log(😄)
,BaseLog=(Base,Num)=>Math.log(Num)/Math.log(Base) // BaseLog_Base(Number)

/*--------------Iterative Taylor/Maclaurin Series--------------*/
/* Note: The functions below all have an equivilant Math API function built into JS, but I consider these as Redundant Syntactic Sugar, so here's what they would look like if we only had the JS syntax we acutally NEEDED for basic functionality

x = number , i = # iterations (more takes more time, but increases accuracy) 

*/

,exp=(x,i)=>{// Exponential: Math.exp(x) | Math.E**x | e^x
    //e^x = 1 + x + x^2/2! + x^3/3! + x^4/4! + ..
    let j=1,N=r=x
    if(i==undefined){// if no i is entered
        let p
        do{p=r;++j;N*=x/j;r+=N}while(p!=r&&![0,NaN,Infinity].includes(r))
        console.log(j-2,'iterations')
    }else if(!Number.isInteger(i)||i<0)return'i must be either zero or a positive integer!'
    else// if i is inputted as a positive integer
        while(j++<=i){N*=x/j;r+=N}
    return Number((r+1).toFixed(14))
}

// Math.E = exp(1) | 2.718281828459045
,EXPM1=(x,i)=>exp(x,i)-1 /* Math.expm1(x) , e.g. EXPM1(1) ➜ 1.71828182845905 , EXPM1(-0.5) ➜ -0.39346934028737
    note that this could've been achieved slightly more efficiently by simply not adding +1 to r, see exp functions above */

/*--------------Other Redundant Syntactic Sugar--------------*/
//Following the Iterative Taylor/Maclaurin Series, here are other simpler Redundant Syntactic Sugar which could've been easily defined in a customised way as opposed to being included in the official syntax:

,IsInt=N=>(N|0)==N /* Number.isInteger(N) , Alternatively: N%1==0, but bitwise operators are faster.
    e.g. IsInt(-3) ➜ true , IsInt(1.2) ➜ false , IsInt(Infinity) ➜ false */

,IsNumeric=N=>!isNaN(N-0) // Basically custom isFinite(N), but IsNumeric(Infinity) ➜ true
,IsFinite=N=>IsNumeric(N)&&N!=Infinity&&N!=-Infinity /* isFinite(N), e.g. IsFinite(Infinity) ➜ false
   IsFinite(-Infinity) ➜ false , IsFinite(NaN) ➜ false , IsFinite(undefined) ➜ false
   IsFinite(Number.MAX_VALUE) ➜ true , IsFinite('0') ➜ true , IsFinite(null) ➜ true
*/
,Sign=N=>N>0?1:N<0?-1:N-0 /* Math.sign(N), e.g. Sign(-0.1) OR Sign(-Infinity) ➜ -1 , Sign(23.4) ➜ OR Sign(Infinity) ➜ 1
    Sign(0) OR Sign(null) ➜ 0 , Sign() OR Sign(undefined) OR Sign('string') OR Sign(NaN) OR Sign([0,1]) ➜ NaN
    Alternatively: !!N|N>>31 but doesn't work for 0 > N > -1, -Infinity, NaN, undefined, strings and arrays */

,ABS=N=>N<0?-1*N:N-0 // Math.abs(N) or |N|, e.g. ABS(-1.5) ➜ 1.5, ABS(2.3) ➜ 2.3
,ABS_Int=N=>{const n=N>>31;return (N^n)-n} /* |N| but faster and only returns integers, 
    e.g. ABS_Int(-1.5) ➜ 1, ABS_Int(2.3) ➜ 2 */

// Methods that turn a fraction into one of two consecutive integers
,Round=N=>N>0?N+0.5|0:IsInt(2*N)?N|0:N-0.5|0 // Math.round(N), e.g. Round(-1.5) ➜ -1, Round(2.3) ➜ 2

,Floor=N=>(N|0)-(N<(N|0)) /* Math.floor(N), Alternatively: (N|0)-(N>>31&1), but doesn't work for 0 > N > -1
ABS(N)==Infinity?N:N<0?(N|0)-1:N|0, but bitwise is faster
    e.g. Floor(0) ➜ 0 , Floor(-4.5) ➜ -5 , Floor(4.5) ➜ 4 */
,Ceiling=N=>(N|0)+(N>(N|0)) /* Math.ceil(N), Alternatively: (N|0)+(N>>31^1), but doesn't work for 0 > N > -1
ABS(N)==Infinity?N:N<0?N|0:(N|0)+1
    e.g. Ceiling(0) ➜ 0 , Ceiling(-4.5) ➜ -4 , Ceiling(4.5) ➜ 5 */

,Trunc=N=>N|0 /* Math.trunc(N), e.g. Trunc(-1.5) ➜ -1 , Trunc(2.3) ➜ 2 , Trunc(-0.1) ➜ 0
    Each of the 7 bitwise operators has a way to truncate on its own to Int32:
    N^0 , ~~N , N&-1 , N<<0 , N>>0 . N>>>0 only works for N >= 0 */
,TruncOut=N=>(N|0)+Sign(N) /* Custom function. Whereas Trunc gets closer to 0, TruncOut gets out and away from 0. 
   e.g. TruncOut(-1.5) ➜ -2 , TruncOut(2.3) ➜ 3  */

,POW=(a,b)=>a**b // Math.pow(a,b), e.g. POW(-2,3) ➜ -8
,IMUL=(a,b)=>(a|0)*(b|0) // Math.imul(a,b), e.g. IMUL(-2.3,3.9) ➜ -6 , basically multiplys the integers

,Max=(a,b)=>a>b?a:b // Math.max(a,b), e.g. Max(0,-Infinity) ➜ 0
,Min=(a,b)=>a<b?a:b // Math.min(a,b), e.g. Min(0,-Infinity) ➜ -Infinity
,Max_Int=(a,b)=>{const n=(a-b)>>31;return b&n|a&~n}/* Max(a,b) but faster and only returns integers, 
    e.g. Max_Int(6.5,-3.5) ➜ 6 */
,Min_Int=(a,b)=>{const n=(a-b)>>31;return a&n|b&~n} /* Min(a,b) but faster and only returns integers, 
    e.g. Min_Int(6.5,-3.5) ➜ -3 */

,SQRT=N=>N**0.5 // Math.sqrt(N) , e.g. SQRT(64) ➜ 8 , Note this result is ±
// Math.SQRT2 = SQRT(2) | 1.4142135623730951 , // Math.SQRT1_2 = SQRT(0.5) | 0.7071067811865476
,CBRT=N=>N**(1/3) // Math.cbrt(N) , e.g. CBRT(64) ➜ 4
,Hypotenuse=(O,A)=>SQRT(O**2+A**2) // Math.hypot(O,A) where O = Opposite & A = Adjacent, e.g. Hypotenuse(3,4) ➜ 5

// Note: you can add 'isNaN(N)?NaN:' at the beginning of each of these, but it's not necessary

,MachinsPI=i=>{// Math.PI = MachinsPI() | 3.141592653589793
    if(i==undefined)i=1000//default iterations
    else if(!Number.isInteger(i)||i<0)return'i must be either zero or a positive integer!'
    let π=0,p,k=-1
    while(p!=π&&++k<i){p=π;π+=(4/(8*k+1)-2/(8*k+4)-1/(8*k+5)-1/(8*k+6))/(16**k)}
    console.log(i,'iterations')
    return π
}

/*--------------Other Operators, Not Native Redundant Syntactic Sugar Methods In JS--------------*/
// Note: anytime a function's name starts with Is, it means it returns true or false

,IsSameSign=(a,b)=>(a^b)>=0 // e.g. SameSign(-43,-5.6) ➜ true , SameSign(0,-5.6) ➜ false

,IsOdd=N=>IsInt(N)&&(N&1)==1 /* true means odd AND indivisible by 2
   e.g. IsOdd(32) ➜ false , IsOdd(33) ➜ true , IsOdd(36) ➜ false , IsOdd(Infinity) ➜ false */
,IsEven=N=>IsInt(N)&&(N&1)==0 /* true means even AND divisible by 2
   e.g. IsEven(2.3) ➜ false , IsEven(0) ➜ true , IsEven(-Infinity) ➜ false */

,IsPowerOf2=N=>IsInt(N)&&N>0&&(N&(N-1))==0 /* e.g. IsPowerOf2(32) ➜ true , IsPowerOf2(36) ➜ false
    Note that techically you can use BaseLog(2,N) and it'll return the exact power, but this is faster.
    BaseLog(B,N) can be used not just for base 2, but also for other bases like decimals */

,IsDivisibleBy=(N,B)=>N%B==0 /* Is N divisible by B? , e.g. IsDivisibleBy(32,2) ➜ true
    IsDivisibleBy(36,2) ➜ true , IsDivisibleBy(33,2) ➜ false , IsDivisibleBy(33,3) ➜ true
    IsDivisibleBy(25,5) ➜ true , IsDivisibleBy(26,5) ➜ false , IsDivisibleBy(25.1,5.2) ➜ false */

,QuotientDivision=(N,D)=>N/D|0 // N = Numerator , D = Denominator, e.g. QuotientDivision(-7,2) ➜ -3

/* Swap the values of two variables such that a=b & b=a
    It's easier to swap without making function, because we're trying to change the original variables, not return new ones

    Temporary Variable: const c=a;a=b;b=c
    Destructuring Assignment: [a,b]=[b,a]
    // these work for every type: fractions|decimal numbers, arrays, strings, boolean, etc
    e.g. let a = '6.5' , b = [-3.5] ; c=a;a=b;b=c OR [a,b]=[b,a] ➜ [ a ➜ [-3.5] ; b ➜ '6.5' ]

    XOR: a^=b;b^=a;a^=b // this only works on integers (or rather only returns integers) but it's faster
    e.g. let a = '6.5' , b = [-3.5] ; a^=b;b^=a;a^=b ➜ [ a ➜ -3 ; b ➜ 6 ]
*/

/*--------------Bitwise and Logical equivalents of each other--------------*/
/* The reason for these functions is that there are 4 logic gate bitwise operators (&, |, ~, ^) 
    and only 3 logical operators (&&, ||, !), so I'm filling the gap here. 
    Note that logical functions return true|false while bitwise functions return 1|0 
*/
,Logical_XOR=(a,b)=>a!=b /* a^b , returns true if the two inputs are unequal, otherwise returns false
    Logical_XOR('Try string','Try string') ➜ false , Logical_XOR('Try string','Another string') ➜ true */

,Bitwise_IsEqual=(a,b)=>2+~(a^b) /* a==b , returns 1 if the two inputs are equal, otherwise returns another integer
    Bitwise_IsEqual(-23.4,8.7) ➜ 32 , Bitwise_IsEqual(-0.5,0) ➜ 1 because both are truncated */

/*--------------Strings: Redundant Syntactic Sugar Methods Native In JS--------------*/
// S = String , C = Character
,TrimStart=(S,C)=>{
    if(C==undefined)C=' ' // if no C is entered, then S.trimStart()
    let s=0;while(S[s]==C)++s
    return S.slice(s)
}// e.g. TrimStart('  Try  ') ➜ 'Try  ' , TrimStart('00034E50.3C0120000','0') ➜ '34E50.3C0120000'
,TrimEnd=(S,C)=>{
    if(C==undefined)C=' ' // if no C is entered, then S.trimEnd()
    let s=-1;while(S.at(s)==C)--s
    return s==-1?S:S.slice(0,s+1)
}// e.g. TrimEnd('  Try  ') ➜ '  Try' , TrimEnd('00034E50.3C0120000','0') ➜ '00034E50.3C012'
,Trim=(S,C)=>{
    if(C==undefined)return TrimStart(TrimEnd(S)) // if no C is entered, then S.trim()
    return TrimStart(TrimEnd(S,C),C)
}/* e.g. Trim('  Try  ') ➜ 'Try' , Trim('00034E50.3C0120000','0') ➜ '34E50.3C012'
    Efficient practices:
        - don't measure .length of whole string if you don't have to
        - use a variable s and slice only once, as opposed to slice once for every character that needs trimming
        - use s+1 as opposed to ++s when you're done using s
        - for Trim, TrimEnd first and only then TrimStart to save on re-indexing characters
*/

,Replace=(SA,P,R)=>{/* Equivilant to the SA.replace(P,R) method but also works with arrays, not just strings.
    Replaces the 1st occurrence of Pattern P with Replacement R. SA is either String '' OR Array []. */
    const i=SA.indexOf(P)
    if(Array.isArray(SA)){//Array
        i==-1||(SA[i]=R);return SA
    }//else if String
        return i==-1?SA:`${SA.slice(0,i)}${R}${SA.slice(i+P.length)}`
}/*e.g. Replace(['S','A','$','S','A','h'],'A','g') ➜ ['S','g','$','S','A','h']
    Replace('SA$SAh','SA','g') ➜ 'g$SAh' , Replace('SA$SAh','','g') ➜ 'gSA$SAh'
    Replace('SA$SAh','D','g') ➜ 'SA$SAh' */
,ReplaceAll=(SA,P,R)=>{/* Equivilant to the SA.replaceAll(P,R) method but also works with arrays, not just strings.
    Replaces ANY & ALL occurrences of Pattern P with Replacement R. SA is either String '' OR Array []. */
    if(Array.isArray(SA)){//Array
        for(let i=SA.indexOf(P);i!=-1;i=SA.indexOf(P,i+1))SA[i]=R
        return SA
    }//else if String
        const p=P.length;if(p==0)return R+SA.split('').join(R)+R
        let i=SA.indexOf(P),AC='' // AC = Already Checked|Changed
        while(i!=-1){
            AC+=SA.slice(0,i)+R;SA=SA.slice(i+p);i=SA.indexOf(P)
        }
        return AC+SA
}/* e.g. ReplaceAll(['S','A','$','S','A','h'],'A','g') ➜ ['S','g','$','S','g','h']
    ReplaceAll('SA$SAh','SA','g') ➜ 'g$gh' , ReplaceAll('SA$SAh','D','g') ➜ 'SA$SAh'
    ReplaceAll('SA$SAh','','g') ➜ 'gSgAg$gSgAghg'
*/

,ToUpperCase=S=>{// S.toUpperCase()
    const L=S.length;let s=''
    for(let i=-1;++i<L;){
        const c=S[i] // c = character
        switch(c){
            case'e':s+='E';break;case't':s+='T';break;case'a':s+='A';break;case'o':s+='O';break
            case'i':s+='I';break;case'n':s+='N';break;case's':s+='S';break;case'h':s+='H';break
            case'r':s+='R';break;case'd':s+='D';break;case'l':s+='L';break;case'c':s+='C';break
            case'u':s+='U';break;case'm':s+='M';break;case'w':s+='W';break;case'f':s+='F';break
            case'g':s+='G';break;case'y':s+='Y';break;case'p':s+='P';break;case'b':s+='B';break
            case'v':s+='V';break;case'k':s+='K';break;case'j':s+='J';break;case'x':s+='X';break
            case'q':s+='Q';break;case'z':s+='Z';break;default:s+=c 
        }
    }
    return s
}// e.g. ToUpperCase('Hello World! אב 🙂 ,./"`@#$%^&*-=+') ➜ 'HELLO WORLD! אב 🙂 ,./"`@#$%^&*-=+'
,ToLowerCase=S=>{// S.toLowerCase()
    const L=S.length;let s=''
    for(let i=-1;++i<L;){
        const c=S[i] // c = character
        switch(c){
            case'E':s+='e';break;case'T':s+='t';break;case'A':s+='a';break;case'O':s+='o';break
            case'I':s+='i';break;case'N':s+='n';break;case'S':s+='s';break;case'H':s+='h';break
            case'R':s+='r';break;case'D':s+='d';break;case'L':s+='l';break;case'C':s+='c';break
            case'U':s+='u';break;case'M':s+='m';break;case'W':s+='w';break;case'F':s+='f';break
            case'G':s+='g';break;case'Y':s+='y';break;case'P':s+='p';break;case'B':s+='b';break
            case'V':s+='v';break;case'K':s+='k';break;case'J':s+='j';break;case'X':s+='x';break
            case'Q':s+='q';break;case'Z':s+='z';break;default:s+=c 
        }
    }
    return s
}// e.g. ToLowerCase('Hello World! אב 🙂 ,./"`@#$%^&*-=+') ➜ 'hello world! אב 🙂 ,./"`@#$%^&*-=+'
/*
Note: there are 2 efficiencies of customising ToUpperCase and ToLowerCase, 
instead of using the Redundant Syntactic Sugar Methods Native In JS of .toUpperCase() & .toLowerCase()

1) you can exclude from this switch letters of other alphabets such as Greek i.e. 'αβδΠΣΩ'
    if you weren't going to use those other non-English letters anyway

2) Sorting out letters in the switch from the most frequently used in literature to the least common,
    as opposed to sorting them alphabetically
*/

/*--------------Strings: Other, Not Redundant Syntactic Sugar Methods Native In JS--------------*/
// Key: S = String

,ReverseCases=S=>{// Replaces capital letters with lowercase and vice versa. Leaves all other characters unchanged.
    const L=S.length;let s='' // s = new string
    for(let i=-1;++i<L;){
        const c=S[i],U=c.toUpperCase() // c = character
        s+=c!=U?U:c.toLowerCase() /* it's more efficient to not even try to convert toUpperCase if it's not needed.
           Capitals were selected first because it is more likely for a character to already be lowercase */
    }
    return s
}//e.g. ReverseCases('Hello World! ,./"`@#$%^&*-=+ αβΣאב❤') ➜ 'hELLO wORLD! ,./"`@#$%^&*-=+ ΑΒσאב❤'

,Acronym=S=>{// Returns the 1st letter of each word in a string S
    if(typeof S!='string')return'Please enter a string'
    S=S.trim()
    let i=S.indexOf(' '),C=S[0] // C = already Checked/Changed

    while(i!=-1){S=S.slice(i+1).trimStart();C+=S[0];i=S.indexOf(' ')}
    return C
}//e.g. Acronym('As Soon   As Possible') ➜ 'ASAP' , Acronym(' Do It Yourself ') ➜ 'DIY'

,NumberOfCharacters=S=>S.length // e.g. NumberOfCharacters('e3ΠΔλθτεαδφξμβζΞΛΨΩΣאב❤') ➜ 23
,ByteSize=S=>new TextEncoder().encode(S).length /* e.g. ByteSize('e3ΠΔλθτεαδφξμβζΞΛΨΩΣאב❤') ➜ 45
Alternatively:
•  ByteSize=S=>new Blob([S]).size , but this requires const {Blob}=require('buffer') ahead of this function if done in NodeJS
•  ByteSize=S=>Buffer.from(S).length , but this only works in NodeJS, not HTML
Both these alternatives do not predict the result ("eager evaluation") in terminal/console, hence TextEncoder().encode() is always preferred
*/
