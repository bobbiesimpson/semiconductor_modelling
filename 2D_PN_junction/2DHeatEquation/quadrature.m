function [ gwt, gpt ] = quadrature( integrationOrder, dimension )

[point1, weight1]=gaussQuad(integrationOrder);

gpt=zeros(integrationOrder^dimension,dimension);
gwt=zeros(integrationOrder^dimension,dimension);

for pt=1:integrationOrder
    index=pt*integrationOrder;
    gpt(index:index+integrationOrder-1,1)=point1(:);
end

