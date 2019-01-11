function [ fbest,pg ] = PSO( population,generation,dimension,low,high)
%particle swarm optimization
%population: #particles, generation: #iterations, dimension 
%low, high: lower and upper bound of variables
%set parameters
w=1;
c1=2;
c2=2;
%pre-allocation of memory
v=zeros(dimension,population);
x=zeros(dimension,population);
pg=zeros(dimension,1);
fibest=-ones(1,population)*inf;
fbest=-inf;
%initialize variable vectors
for i=1:population
    x(:,i)=rand(1,dimension).*(high-low)+low;
end
pi=x;
for i=1:generation
    for j=1:population
        f=calculate_fitness(x(:,j));%calculate objective function
        if(f>fibest(j))
            fibest(j)=f;
            pi(:,j)=x(:,j);
        end
        if(f>fbest)
            fbest=f;%update globla optimum
            pg=x(:,j);%update local optimum
        end
    end
    %update the speed and location of particles
    for j=1:population
        v(:,j)=w*v(:,j)+c1*rand(dimension,1).*(pi(:,j)-x(:,j))+c2*rand(dimension,1).*(pg-x(:,j));%update speed of particles
        x(:,j)=x(:,j)+v(:,j);%update the location of particles
    end
end
end

