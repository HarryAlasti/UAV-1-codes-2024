function [Pos , dist, firstpnt, num_pnts] = Parser3(CTR)

Pos = [];
firstpnt = [];
num_pnts = [];
dist = [];

C = CTR';
[row , col] = size(C);
k = row;
blk = 1;
while k >= 1
    cnt = C(1 , 2);
    num_pnts = [num_pnts, C(1,2)];
    firstpnt = [firstpnt ; [C(2, 1)  C(2,2)   C(1,1)]];
    Block = C(2 : cnt , :);
    CBlock = [Block(1:end , :) , C(1,1) * ones(cnt-1 , 1) , blk * ones(cnt-1, 1)];
    
    dist = [dist , sum(sqrt((Block(1:end-1,1) - Block(2:end,1)).^2 + (Block(1:end-1,2) - Block(2:end,2)).^2))];
    
    Pos = [Pos ; CBlock];
    C = C(cnt + 2 : end , :);
    k = k - (cnt + 1);
    blk = blk + 1;
end
return,