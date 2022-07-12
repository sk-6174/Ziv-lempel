clc
clear
close all
%19EC30030
fileID = fopen("MnM.txt","r");
A = fscanf(fileID,"%c",Inf);
A(A==' ') =[];
A = erase(A,newline);
len_A = length(A);
disp("Original Text")
disp(A) %original tex
%encoding for 20622 characters which does not include spaces or newline
% '•' before second occurance of Our Brand(line 463 in text) is conidered distinct from '.'
K = ['!', '"',"'", ',', '-', '.', '4', ':', ';', '?', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'R', 'S', 'T', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '•'];
M = length(K); 
bits = ceil(log2(M)); %no. of bits for each character
V = 0:1:M-1;
V = dec2bin(V); %59x6 

A_map = containers.Map;

for i=1:M
A_map(K(i)) = V(i,:);
end

A_map2 = containers.Map;
for i=1:M
A_map2(V(i,:)) = K(i);
end
var = A_map2;

W=2000;
Pi = 1;
Pf = W;

initial = '';
for i = 1:W
    initial = strcat(initial,A_map(A(i)));
end

while Pf < len_A
    [encoded, n] = encode(Pi,Pf,A,W,A_map);
    Pi = Pi + n;
    Pf = Pf + n;
    initial = strcat(initial, encoded);
end
encoded_string = initial; %Encoded
disp("Encoded:")
disp(encoded_string)
init_decoded = decode1(encoded_string, W, var);
initial = strcat(init_decoded,initial(6*W+1:length(initial)));

Pointer = 1;
st = '';
req_st = encoded_string((6*W+1):length(encoded_string));

while Pointer < length(req_st)
    [decoded_, n] = decode(req_st, W, var, Pointer, init_decoded);
    Pointer = Pointer + n;
    decoded_ = char(decoded_);
    init_decoded = strcat(init_decoded,decoded_);
end
decoded_text = init_decoded; % Decoded
disp("Decoded:")
disp(decoded_text)

Commpression_ratio = (len_A*bits-length(encoded_string))/(len_A*bits);
disp("Commpression ratio:")
disp(Commpression_ratio)





function [code, count]= encode(Pii, Pfi, A, Win_size,A_map)
    stop = 1;
    count = 0;
    
    while stop
        if Pfi+count+1 > length(A)
            break
        else
        stop = contains(A(Pii:Pfi), A((Pfi+1):(Pfi+count+1)));
        count = count + 1;
        end
    end
    count=count-1;
    if count > 1
        match = A((Pfi+1):(Pfi+count));
        ind = strfind(A(Pii:Pfi),match);
        ind = ind(length(ind));
        ind = Win_size-ind+1;
        position = ind;
        num =  count;
        num = dec2bin(num);
        sn = 2*length(num)-1;
        num = count;
        num = dec2bin(num, sn);
        un = ceil(log2(Win_size));
        num2 = position;
        num2 = dec2bin(num2,un);
        code = strcat(num,num2);

    elseif count == 0
        match = A(Pfi+1);
        neau = A_map(match);
        code = strcat('1',neau);
        count = count+1;
    elseif count == 1
        match = A(Pfi+1);
        neau = A_map(match);
        code = strcat('1',neau);
    else
        code = '';
    end
end


function decoded1 = decode1(coded, Win, new)
    coded_in = coded(1:(6*Win));
    decoded1 = '';
    for i = 1:Win
        code_sample = coded_in(1+6*(i-1):6*i);
        ch = new(code_sample);
        decoded1 = strcat(decoded1, ch);
        decoded1 = char(decoded1);
    end
end

function [decoded, pos] = decode(coded, Win, new,P, decoded_init)
    code_in = coded; % req_st
    u_size = ceil(log2(Win));
    l = length(decoded_init);

    %decoded_init is decoded string till then
    if code_in(P:P+4) == '00000'
        k = 5;
        nin = code_in((P+k):(P+(2*k)));
        nin = bin2dec(nin);
        uin = code_in((P+1+(2*k)):(P+(2*k)+u_size));
        uin = bin2dec(uin);
        match = decoded_init((l-(uin-1)):(l-uin+nin));
        decoded = match;
        pos = (2*k) + u_size +1;


    elseif code_in(P:P+3) == '0000'
        k = 4;
        nin = code_in((P+k):(P+(2*k)));
        nin = bin2dec(nin);
        uin = code_in((P+1+(2*k)):(P+(2*k)+u_size));
        uin = bin2dec(uin);
        match = decoded_init((l-(uin-1)):(l-uin+nin));
        decoded = match;
        pos = (2*k) + u_size +1;
    
    elseif code_in(P:P+2) == '000'
        k = 3;
        nin = code_in((P+k):(P+(2*k)));
        nin = bin2dec(nin);
        uin = code_in((P+1+(2*k)):(P+(2*k)+u_size));
        uin = bin2dec(uin);
        match = decoded_init((l-(uin-1)):(l-uin+nin));
        decoded = match;
        
        pos = (2*k) + u_size +1;
     
    elseif code_in(P:P+1) == '00'
        k = 2;
        nin = code_in((P+k):(P+(2*k)));
        nin = bin2dec(nin);
        uin = code_in((P+1+(2*k)):(P+(2*k)+u_size));
        uin = bin2dec(uin);
        match = decoded_init((l-(uin-1)):(l-uin+nin));
        decoded = match;
        pos = (2*k) + u_size +1;

    
    elseif code_in(P) == '0'
        k = 1;
        nin = code_in((P+k):(P+(2*k)));
        nin = bin2dec(nin);
        uin = code_in((P+1+(2*k)):(P+(2*k)+u_size));
        uin = bin2dec(uin);
        match = decoded_init((l-(uin-1)):(l-uin+nin));
        decoded = match;
        pos = (2*k) + u_size +1;
     
    elseif code_in(P) == '1'
        uin = code_in(P+1:P+6);
        ch = new(uin);
        pos = 7;
        decoded = ch;
    end
    
end

