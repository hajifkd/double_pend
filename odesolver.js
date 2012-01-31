function add_vector(a,x,b,y){
    var result=new Array(x.length);
    for(var i=0;i<x.length;i++){
        result[i]=a*x[i]+b*y[i];
    }
    return result;
}

function runge_kutta2(dt,vec,func){
    var eta=add_vector(1,vec,dt/2,func(vec));
    return add_vector(1,vec,dt,func(eta));
}

function euler(dt,vec,func){
    return add_vector(1,vec,dt,func(vec));
}

function symplectic1(dt, vec, func){
    var tmp = new Array(vec.length);
    var res1 = func(vec);
    var i;


    for(i = 0; i < vec.length / 2; i++){
        tmp[i] = vec[i] + res1[i] * dt;
    }

    for(; i< vec.length; i++){
        tmp[i] = vec[i];
    }

    var res2 = func(tmp);
    var result = new Array(vec.length);

    for(i = 0; i < vec.length / 2; i++){
        result[i] = tmp[i];
    }

    for(; i< vec.length; i++){
        result[i] = tmp[i] + res2[i] * dt;
    }

    return result;
}

function runge_kutta4(dt,vec,func){
    var k1=func(vec);
    var k2=func(add_vector(1,vec,dt/2,k1));
    var k3=func(add_vector(1,vec,dt/2,k2));
    var k4=func(add_vector(1,vec,dt,k3));
    var tmp=new Array(vec.length);
    for(var i=0;i<vec.length;i++){
        tmp[i]=vec[i]+dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }
    return tmp;
}

var oldx;

function init_leap_frog2(dt,vec,func){
    oldx=add_vector(1,vec,-dt/2,func(vec));
}

function __leap_frog2(dt,vec,func){
    var nv=add_vector(1,oldx,dt,func(vec));
    oldx=vec;
    return nv;
}

function leap_frog2(dt,vec,func){
    return __leap_frog2(dt,__leap_frog2(dt,vec,func),func);
}

var old2;
var old1;

function init_adams_bashforth3(dt,vec,func){
    old1=func(runge_kutta4(-dt,vec,func));
    old2=func(runge_kutta4(-dt,old1,func));
}

function adams_bashforth3(dt,vec,func){
    var old0=func(vec);
    var tmp=new Array(vec.length);
    for(var i=0;i<vec.length;i++){
        tmp[i]=vec[i]+dt/12*(23*old0[i]-16*old1[i]+5*old2[i]);
    }
    old2=old1;
    old1=old0;
    return tmp;
}

var ad_old=new Array(5);

function init_adams_bashforth6(dt, vec, func){
    ad_old[0] = func(runge_kutta4(-dt,vec,func));
    for(var i=1;i<ad_old.length;i++){
        ad_old[i] = func(runge_kutta4(-dt, ad_old[i-1], func));
    }
}

function adams_bashforth6(dt,vec,func){
    var old0=func(vec);
    var tmp=new Array(vec.length);
    for(var i=0;i<vec.length;i++){
        tmp[i]=vec[i]+dt*(4277*old0[i] - 7923*ad_old[0][i] + 9982*ad_old[1][i]
                - 7298*ad_old[2][i] + 2877*ad_old[3][i] - 475*ad_old[4])/1440;
    }
    for(var i=1;i<ad_old.length;i++){
        ad_old[i] = ad_old[i-1];
    }
    ad_old[0]=old0;
    return tmp;
}
