function [ GlobalBest, BestCost ] = MEPSO( param , model ,fobj  )

%PSO Infotmation

c=0;
d=1;
dim=model.nVar; 
SolbestX=struct();   
Vmax=6;
noP=param.nPop;
wMax=0.6;
wMin=0.2;
c1=2;
c2=2;
pop=noP;
% Initializations
iter=param.MaxIt;
vel=zeros(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
% cg_curve=zeros(1,iter);
lb= c.*ones( 1,dim );    % 下界
ub= d.*ones( 1,dim );    % 上届
% Random initialization for agents.
pos=initialization(noP,dim,ub,lb); 


for i=1:noP
    pBestScore(i)=inf;
end

% Initialize gBestScore for a minimization problem
 gBestScore=inf;
     
% for i = 1 : noP
    for i=1:size(pos,1)   
        Flag4ub=pos(i,:)>ub;
     Flag4lb=pos(i,:)<lb;
     pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %Calculate objective function for each particle
        [fitness( i ),Sol(i)] = fobj( pos( i, : ) ) ;

        if(pBestScore(i)>fitness(i))
            pBestScore(i)=fitness(i);
            pBest(i,:)=pos(i,:);
        end
        if(gBestScore>fitness(i))
            gBestScore=fitness(i);
            gBest=pos(i,:);
            SolbestX=Sol(i);
            Row_Id=i;
        end
     end                       
% end  

for l=1:iter 
    threshold=2-l*((2)/iter);
  
    w=wMax-l*((wMax-wMin)/iter);
    %Update the Velocity and Position of particles
    for i=1:size(pos,1)
        for j=1:size(pos,2)       
            vel(i,j)=w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+c2*rand()*(gBest(j)-pos(i,j));
            
            if(vel(i,j)>Vmax)
                vel(i,j)=Vmax;
            end
            if(vel(i,j)<-Vmax)
                vel(i,j)=-Vmax;
            end            
            pos(i,j)=pos(i,j)+vel(i,j);
        end
    end
   
   
     for i=1:size(pos,1)   
        Flag4ub=pos(i,:)>ub;
     Flag4lb=pos(i,:)<lb;
     pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %Calculate objective function for each particle
        [fitness( i ),Sol(i)] = fobj( pos( i, : ) ) ;

        if(pBestScore(i)>fitness(i))
            pBestScore(i)=fitness(i);
            pBest(i,:)=pos(i,:);
        end
        if(gBestScore>fitness(i))
            gBestScore=fitness(i);
            gBest=pos(i,:);
            SolbestX=Sol(i);
            Row_Id=i;
        end
     end

     for j=1:dim
          LB(:,j)=min(pos(:,j))  ;
          UB(:,j)=max(pos(:,j));
     end
 if mod(iter,50)==0      
    
     m=size(pos,1); n=dim;
a=2;
k=-1/a;
I=k*log(1-rand(m, n));%指数分布噪声exponential distribution noise
qq=fitnessDistanceBalance( pos, fitness );
 gg=pos(qq,:);%Individuals obtained from FDB
 ak=2*rand(pop, dim) - 1;%random matrix between [-1,1] Q
 for i=1:size(pos,1)  
pop_x(i,:)=0.5*(UB-LB).*I(i,:)+ak(i,:).*gg;

         
        Flag4ub=pop_x(i,:)>ub;
     Flag4lb=pop_x(i,:)<lb;
     pop_x(i,:)=(pop_x(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        [fitness1( i ),Sol(i)] = fobj( pop_x( i, : ) ) ;
        if(gBestScore>fitness1(i))
            gBestScore=fitness1(i);
            gBest=pop_x(i,:);
            SolbestX=Sol(i);
            Row_Id=i;
        end
        if fitness1( i )<fitness(i)
            
            pos(i,:)=pop_x(i,:);
        end
 end
    
 else   

   [pos1]=corOppose2(pos,UB,LB,threshold,Row_Id);

 for i=1:size(pos,1)   
        Flag4ub=pos1(i,:)>ub;
     Flag4lb=pos1(i,:)<lb;
     pos1(i,:)=(pos1(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %Calculate objective function for each particle
        [fitness2( i ),Sol(i)] = fobj( pos1( i, : ) ) ;
        if(gBestScore>fitness(i))
            gBestScore=fitness(i);
            gBest=pos1(i,:);
            SolbestX=Sol(i);
            Row_Id=i;
        end
         if fitness2( i )<fitness(i)
            
            pos(i,:)=pos1(i,:);
          end
 end 
 end

    BestCost(l)=gBestScore;
    GlobalBest.Sol=SolbestX;
    disp(['Iteration ' num2str(l) ': Best Cost = ' num2str(BestCost(l))]);
    GlobalBest.BestCost = BestCost(1:l);
    GlobalBest.MaxIt =iter; 
end

end
function [Positions]=corOppose2(Positions,upper,lower,threshold,Row_Id)

[n,b] = size(Positions);
for i=1:n(1) 
    
    if i ~= Row_Id
        
        sum=0;
        greater=[];
        less=[];
        x=1;z=1;y=1;
        
        for j=1:b
            d(x)=abs(Positions(Row_Id,j)-Positions(i,j));
            if d(x)<threshold
                greater(y)=j;
                y=y+1;
            else
                less(z)=j;
                z=z+1;
            end
            sum=sum+d(x)*d(x);
            x=x+1;
        end
        %     sum
        src=1-(double(6*sum))/(double(n(1)*(n(1)*n(1)-1)));
        %     src
        if src<=0
            if size(greater)<size(less)
                %             for j=1:size(less)
                %                 dim=less(j);
                %                 Positions(i,dim)=ub(dim)+lb(dim)-Positions(i,dim);
                %             end
            else
                for j=1:size(greater)
                    dim=greater(j);
                    Positions(i,dim)=(upper(1,dim)+lower(1,dim)-Positions(i,dim));
                end
            end
        end
    end
end
end