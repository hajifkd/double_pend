var epsilon=1e-8

function partial_diff3(f,x,i){
    var x1=x.slice();
    var x2=x.slice();
    x1[i]+=epsilon;
    x2[i]-=epsilon;
    return (f(x1)-f(x2))/(2*epsilon);
}

function partial_diff5(f,x,i){
    var x1=x.slice();
    var x2=x.slice();
    var x3=x.slice();
    var x4=x.slice();
    x1[i]-=2*epsilon;
    x2[i]-=epsilon;
    x3[i]+=epsilon;
    x4[i]+=2*epsilon;
    return (f(x1)-8*f(x2)+8*f(x3)-f(x4))/(12*epsilon);
}

var partial_diff=partial_diff5;

function gen_canonical_eq(h){
    return function(vec){
        var result=new Array(vec.length);
        var len=vec.length/2;
        for(var i=0;i<len;i++){
            result[i]=partial_diff(h,vec,len+i);
            result[i+len]=-partial_diff(h,vec,i);
        }
        return result;
    };
}
