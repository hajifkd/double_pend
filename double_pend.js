var svg;
var canvas;

var line1;
var line2;

var circle1;
var circle2;

var fpos1;
var fpos2;

var pos1;
var pos2;

var vec;
var ratio;

var timer;

var solver=runge_kutta4;//odeのソルバー。runge_kutta{2,4},leap_frog2,euler

var centerx=400;//中心のx
var centery=300;//中心のy
var radius1=120;//動径長(cm)
var l1=radius1*0.01;

var radius2=120;//動径長(cm)
var l2=radius2*0.01;

var m1=0.1;//質量(kg)
var m2=0.1;//質量(kg)

var g=9.80665;//重力加速度(m/s^2)

var t=0;
var dt=0.0001;

var draw=true;

var H=function(x){
    return (l2*l2*m2*x[2]*x[2]
            +l1*l1*(m1+m2)*x[3]*x[3]
            -2*m2*l1*l2*x[2]*x[3]*Math.cos(x[0]-x[1]))/
        (2*l1*l1*l2*l2*m2*(m1+Math.pow(Math.sin(x[0]-x[1]),2)*m2))
        -m2*g*l2*Math.cos(x[1])-(m1+m2)*g*l1*Math.cos(x[0]);
};//Hamiltonian H(q,p)

var func=gen_canonical_eq(H);//canonical eq

function rad2pos(radius,radian,stdx,stdy){
    var x=radius*Math.sin(radian);
    var y=radius*Math.cos(radian);
    return [x+stdx,y+stdy];
}

function start(){
    var gr=svg.group({stroke:'blue', strokewidth: 2});
    pos1=rad2pos(radius1,vec[0],centerx,centery);
    pos2=rad2pos(radius2,vec[1],pos1[0],pos1[1]);
    canvas.moveTo(pos2[0],pos2[1]);
    canvas.beginPath();
    timer=setInterval(function(){
        if(draw) canvas.moveTo(pos2[0],pos2[1]);
        for(var j=0;j<20*ratio;j++){
            for(var i=0;i<15;i++){
                t+=dt;
                vec=solver(dt,vec,func);
            }
            fpos1=pos1;
            fpos2=pos2;
            pos1=rad2pos(radius1,vec[0],centerx,centery);
            pos2=rad2pos(radius2,vec[1],pos1[0],pos1[1]);
            if(draw) canvas.lineTo(pos2[0],pos2[1]);
        }
        $("#energy").text(H(vec).toFixed(8));
        $("#time").text(t.toFixed(5));
        if(draw){
            canvas.stroke();
            canvas.beginPath();
        }
        svg.change(line1,{x2:pos1[0],y2:pos1[1]});
        svg.change(line2,{x1:pos1[0],y1:pos1[1],x2:pos2[0],y2:pos2[1]});
        svg.change(circle1,{cx:pos1[0],cy:pos1[1]});
        svg.change(circle2,{cx:pos2[0],cy:pos2[1]});
    },30);
}

function init(){
    canvas.clearRect(-10,-10,810,610);
    var i_q1=parseFloat($("#initq1").val());
    var i_p1=parseFloat($("#initp1").val());
    var i_q2=parseFloat($("#initq2").val());
    var i_p2=parseFloat($("#initp2").val());
    var i_m1=parseFloat($("#m1").val());
    var i_m2=parseFloat($("#m2").val());
    var i_g=parseFloat($("#g").val());
    var i_l1=parseFloat($("#l1").val());
    var i_l2=parseFloat($("#l2").val());
    var i_r=parseFloat($("#ratio").val());

    if(isNaN(i_q1)) i_q1=0;
    if(isNaN(i_p1)) i_p1=0;
    if(isNaN(i_q2)) i_q2=0;
    if(isNaN(i_p2)) i_p2=0;
    vec=[i_q1,i_q2,i_p1,i_p2];

    if(i_m1>0) m1=i_m1;
    if(i_m2>0) m2=i_m2;
    if(i_g>=0) g=i_g;
    if(i_l1>0){
        l1=i_l1;
        radius1=l1*100;
    }
    if(i_l2>0){
        l2=i_l2;
        radius2=l2*100;
    }

    ratio=1;
    if(i_r>0) ratio=i_r;

    draw=$("#draw").attr("checked");

    switch($("#ode").val()){
        case "rk4":
            solver=runge_kutta4;
            break;
        case "rk2":
            solver=runge_kutta2;
            break;
        case "lf":
            solver=leap_frog2;
            init_leap_frog2(dt,vec,func);
            break;
        case "euler":
            solver=euler;
            break;
        case "ab3":
            solver=adams_bashforth3;
            init_adams_bashforth3(dt,vec,func);
            break;
        case "ab6":
            solver=adams_bashforth6;
            init_adams_bashforth6(dt,vec,func);
            break;
        case "symplectic1":
            solver = symplectic1;
            break;
    }

    clearInterval(timer);
    $("#tenergy").text(H(vec).toFixed(8));
    t=0;
    
    pos1=rad2pos(radius1,vec[0],centerx,centery);
    pos2=rad2pos(radius2,vec[1],pos1[0],pos1[1]);
    svg.change(line1,{x2:pos1[0],y2:pos1[1]});
    svg.change(line2,{x1:pos1[0],y1:pos1[1],x2:pos2[0],y2:pos2[1]});
    svg.change(circle1,{cx:pos1[0],cy:pos1[1]});
    svg.change(circle2,{cx:pos2[0],cy:pos2[1]});
}

$(function(){
    $("#svgcanvas").svg();
    svg=$("#svgcanvas").svg('get');
    canvas=$("#canvas")[0].getContext('2d');
    canvas.lineWidth=0.5;
    canvas.strokeStyle="green";
    var gr=svg.group({stroke:'red', strokewidth: 2});
    vec=[2,1,0,0];//(q,p)
    pos1=rad2pos(radius1,vec[0],centerx,centery);
    pos2=rad2pos(radius2,vec[1],pos1[0],pos1[1]);
    line1=svg.line(gr,centerx,centery,pos1[0],pos1[1]);
    line2=svg.line(gr,pos1[0],pos1[1],pos2[0],pos2[1]);
    circle1=svg.circle(pos1[0],pos1[1],10,{fill:'orange'});
    circle2=svg.circle(pos2[0],pos2[1],10,{fill:'skyblue'});
    $("#tenergy").text(H(vec).toFixed(8));
    $("#initq1").val(vec[0]);
    $("#initp1").val(vec[2]);
    $("#initq2").val(vec[1]);
    $("#initp2").val(vec[3]);
    $("#m1").val(m1);
    $("#m2").val(m2);
    $("#g").val(g);
    $("#l1").val(l1);
    $("#l2").val(l2);
    $("#ratio").val(1);
    $("#sb").click(function(){
        init();
        start();
    });

    $("#stopbutton").click(function(){
        clearInterval(timer);
    });

    $("#resumebutton").click(function(){
        start();
    });

    $("#initbutton").click(function(){
        init();
    });

    //start(vec,1);
});
