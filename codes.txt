%Simplex--------------------------------------------------------------------------
clc
clear all
format short
no_of_var=2;
C=[2 5]
a=[1 4; 3 1; 1 1]
b=[24; 21; 9]
s=eye(size(a,1))
A=[a s b] 
cost=zeros(1,size(A,2))
cost(1:no_of_var)=C
bv=no_of_var+1:1:size(A,2)-1
zjcj=cost(bv)*A-cost
zcj=[zjcj;A]

simpletable=array2table(zcj);
simpletable.Properties.VariableNames(1:size(zcj,2))={'x1','x2','s1','s2','s3','sol'}

RUN=true;
while RUN
    if any(zjcj<0);
        fprintf('The current BFS is not optimal')
        zc=zjcj(1:end-1)
        [Enter_var,pvt_col]=min(zc)
        if all(A(:,pvt_col)<=0)
            error('LPP is Unbounded')
        else
            sol=A(:,end)
            column=A(:,pvt_col)

            for i=1:size(A,1)
                if column(i)>0
                    ratio(i)=sol(i)./column(i)
                else
                    ratio(i)=inf
                end
            end
            [leaving_val,pvt_row]=min(ratio)
        end
        bv(pvt_row)=pvt_col
        pvt_key=A(pvt_row,pvt_col)
        A(pvt_row,:)=A(pvt_row,:)/pvt_key
        for i=size(A,1)
            if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row)
            end
        end
        zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:)
        zcj=[zjcj;A]
        table=array2table(zcj)
        table.Properties.VariableNames(1:size(zcj,2))={'x1','x2','s1','s2','s3','sol'}
    else
        RUN=false
        fprintf('The current BFS is Optimal solution')
    end
end



%Dual Simplex-----------------------------------------------------------------------------------------------------

variables = {'x1', 'x2', 'x3', 's1', 's2', 'sol'}
cost = [-2 0 -1 0 0 0]
info = [-1 -1 1; -1 2 -4]
b = [-5; -8]
s = eye(size(info, 1))
A = [info s b]
BV = []

for j = 1 : size(s, 2)
    for i = 1 : size(A, 2)
        if A(:, i) == s(:, j)
            BV = [BV i]
        end
    end
end
fprintf("Basic varibales (BV) = ")
disp(variables(BV))

ZjCj = cost(BV)*A - cost
ZCj = [ZjCj; A]
simpleTable = array2table(ZCj)
simpleTable.Properties.VariableNames(1: size(ZCj, 2)) = variables

RUN =  true
while RUN
    sol = A(:, end)
    if any(sol < 0)
        fprintf('current bfs is not feasible sol')
        [Leaving.value, Pvt_Row] = min(sol)
        fprintf('Leaving Row = %d\n', Pvt_Row)
        Row = A(Pvt_Row, 1:end-1)
        ZRow = ZjCj(:, 1:end-1)
        for i = 1 : size(Row, 2)
            if Row(i) < 0
                ratio(i) = abs(ZRow(i)/Row(i))
            else
                ratio(i) = inf
            end
        end
        [MinRatio, Pvt_Col] = min(ratio)
        fprintf('Entering Variable is %d\n', Pvt_Col)
        BV(Pvt_Row) = Pvt_Col
        fprintf('Basic variables(BV) = ')
        disp(variables(BV))

        Pvt_Key = A(Pvt_Row, Pvt_Col)
        A(Pvt_Row, :) = A(Pvt_Row, :)/ Pvt_Key
        for i = 1 : size(A, 1)
            if i ~= Pvt_Row
                A(i, :) = A(i, :) - A(i, Pvt_Col)*A(Pvt_Row, :)
            end
            ZjCj = ZjCj - ZjCj(Pvt_Col)*A(Pvt_Row, :)
            ZCj = [ZjCj; A];
            simpleTable = array2table(ZCj)
            simpleTable.Properties.VariableNames(1: size(ZCj, 2)) = variables
        end
    else
        RUN = false
        fprintf('current bfs is feasible and optimal\n')
    end
end

% final optimal solution
Final_BFS = zeros(1, size(A, 2))
Final_BFS(BV) = A(:, end)
Final_BFS(end) = sum(Final_BFS.*cost)
optimal_BFS = array2table(Final_BFS)
optimal_BFS.Properties.VariableNames(1: size(optimal_BFS, 2)) = variables



% LCM Method (Least Cost method)-----------------------------------------------------------------------------------

clc
clear all
C = [6 4 1 5; 8 9 2 7; 4 3 6 4]
a = [14 16 7]
b = [6 10 15 4]
m = size(C, 1)
n = size(C, 2)

if sum(a) == sum(b)
    fprintf('given problem is balanced\n');
else
    fprintf('given problem is not balanced\n');
    if sum(a) < sum(b)
        C(end+1, :) = zeros(1, length(b))
        a(end+1) = sum(b) - sum(a)
    else
        C(:, end+1) = zeros(length(a), 1)
        b(end+1) = sum(a) - sum(b)
    end
end

x = zeros(m, n);
Initial_C = C;
Z = 0;

for i = 1 : size(C, 1)
    for j = 1 : size(C, 2)
        Cpq = min(C(:));
        if Cpq == Inf
            break
        end
        [p1, q1] = find(Cpq == C)
        Xpq = min(a(p1), b(q1));
        [x(p1, q1), ind] = max(Xpq)
        P = p1(ind)
        Q = q1(ind)
        X(P, Q) = min(a(P), b(Q))

        if(min(a(P), b(Q)) == a(P))
            b(Q) = b(Q) - a(P)
            a(P) = a(P) - x(P, Q)
            C(P, :) = Inf;
        else
            a(P) = a(P) - b(Q)
            b(Q) = b(Q) - X(P, Q)
            C(:, Q) = Inf;
        end
    end
end

for i = 1 : size(C, 1)
    for j = 1 : size(C, 2)
        Z = Z + Initial_C(i, j)*X(i, j)
    end
end

array2table(X)
fprintf('transporatation cost is %d', Z)


% Stepest Descent Method------------------------------------------------------------------------------------------
% f(x1,x2)= x1-x2+2x1^2+2x1x2+x2^2


clc
clear all
syms x1 x2
% Define objective function (min)
f1= x1-x2+2*x1^2+2*x1*x2+x2^2;
%f1=126*x1+182*x2-9*x1^2-13*x2^2;
% f1=10*(x2-x1^2)+(1-x1)^2
fx=inline(f1)    % convert to function
fobj=@(X) fx(X(1),X(2))

% Gradient of f
grad=gradient(f1)
G1=inline(grad)
Gx=@(X) G1(X(1),X(2))

% Hessian matrix
H1=hessian(f1)
Hx=inline(H1)

X0=[0 0]     % Set initial vector
maxiter=10;   % Set max iteration
tol=10^(-3); % Max tolerance
iter=0;      % Initial counter

% Steps
X=[];
while norm(Gx(X0))>tol && iter<maxiter
X=[X;X0]      % Save all vectors
S=-Gx(X0)     % Compute gradient at X
H=Hx(X0)      % Compute Hessian at X
lambda=S'*S/(S'*H*S)  % Compute lambda
Xnew=X0+lambda*S'     % Update X
X0=Xnew               % Save new X
iter=iter+1           % Update iteration
end

fprintf('Optimal solution is (%f, %f)\n', X0(1),X0(2))
fprintf('Optimal value is %f', fobj(X0))
