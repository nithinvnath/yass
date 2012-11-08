% Decoding a bit stream from a stego image.

img = imread('D:\newimage.jpg');
seed = load('seed.mat','s');
rng(seed.s);

q = 7;           % Rate of enccoding for RA Codes
len = 8;        % Length of encoded string
permvector = randperm(len*8*q ); %Finding the same random permutation used for encoding
%Finding inverse permuation vector
ipermvector = permvector;
for i=1:len*8*q
    ipermvector(permvector(i)) = i;
end

[M N] = size(img);
B = 32;                     %B = block size
Mb = floor(M/B);            
Nb = floor(N/B);            %Mb x Nb B sized blocks


%-------------Finding Quantization Matrix-----------------------------%
Qf = 50;                    %Quality factor for quantization
%Quantization Table
QM=[16  11  10  16  24   40   51   61   
12  12  14  19  26   58   60   55   
14  13  16  24  40   57   69   56   
14  17  22  29  51   87   80   62   
18  22  37  56  68   109  103  77   
24  35  55  64  81   104  113  92   
49  64  78  87  103  121  120  101  
72  92  95  98  112  100  103  99];

%Scaling quantization table based on Qf
if Qf<50
    S = 5000/Qf;
else
    S = 200 - 2*Qf;
end
newQM=QM;
if S~=0
    newQM = floor((QM * S + 50)/100);
end


%--------------------------------------------------------------------------%
bitstream = [];
for I = 1:Mb
    for J = 1:Nb
        %Selecting the random 8x8 block from BxB block
        Sx = randi([1,B-7]);
        Sy = randi([1,B-7]);
        randBlock = double(img((I-1)*B+Sx:(I-1)*B+7+Sx,(J-1)*B+Sy:(J-1)*B+7+Sy));
        %Computing DCT for the random block
        randBlock = randBlock-128;
        randBlockDCT = dct2(randBlock);
        randBlock = randBlock+128;
        randBlockQ = (randBlockDCT ./ newQM);
        colstep = 1; rowstep = -1;
        i=1;j=1;
        %travelling zigzag
        for n = 2:20;       %first 19 AC coefficients
           
            i=i+rowstep;
            j=j+colstep;
            if (i<1)
                i=i+1;
                colstep=colstep*-1;
                rowstep=rowstep*-1;
            end
            if (j<1)
                j=j+1;
                colstep=colstep*-1;
                rowstep=rowstep*-1;
            end
           absval = abs(round(randBlockQ(i,j)));
            if( absval > 0) 
                if(mod(absval,2) == 0)
                    bitstream = [bitstream, 0];
                elseif(mod(absval,2) ~=0)
                   bitstream = [bitstream, 1];
               end
            end
        end
        
    end
end

% Getting the bitstream by reapplying the permutation i.e. applying the
% inverse permutation

bitstream = bitstream(1:len*8*q);
drepbitstream = bitstream;          %Preallocating 
for i=1:len*8*q 
    drepbitstream(1,i) = bitstream(ipermvector(1,i));        %getting back repeated bit stream
end

% Finding the correct bit per block   
obitstream = double(len*8);
bc=1;
for i=1:7:len*8*q
   if(histc(drepbitstream(i:i+6),1) >= 5)
	obitstream(bc) = 1;
    bc=bc+1;
   elseif(histc(drepbitstream(i:i+6),0) >=5)
    obitstream(bc) = 0;
    bc=bc+1;
   else
       i
       drepbitstream(i:i+6)
   end
end


decodestr = char(len);
cc = 1;
for i=1:8:len*8
    decodestr(cc) = char(bin2dec(dec2hex(obitstream(i:i+7))'));
    cc = cc+1;
end

decodestr
    

