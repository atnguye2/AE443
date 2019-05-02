function cost = value(N)
cost = 2*ones(size(N));
cost(N==0)=4;
cost(N==1)=3;
cost(N==2)=2;
cost(N==3)=1;
cost(N==4)=0;
cost(N==5)=.25;
cost(N==6)=.5;
cost(N==7)=.75;
cost(N==8)=1;
end