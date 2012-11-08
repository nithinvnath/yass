img = imread('image.jpg');

seed = load('seed.mat','s');
rng(seed.s);


%------------------------------------------------------------------------
str = 'qwertyui';
strbin = dec2bin(str,8);
bitstream = reshape(strbin',1,length(str)*8);
bitstream = bitstream - '0';

%Repeat Accumulate Encoding

q = 7;      % Rate of enccoding for RA Codes
repbitstream = zeros(q, length(str)*8);
for i=1:q
repbitstream(i,:) = bitstream;
end

repbitstream = repbitstream(:)';      %Repeating Bits.

% Finding a random permutation of bitsream %
permvector = randperm(length(str)*8*q );
permbitstream = repbitstream;   %Preallocating 
for i=1:length(str)*8*q 
    permbitstream(1,i) = repbitstream(permvector(1,i));
end

% Accumulating bitstream
accbitstream = uint8(zeros(1,length(str)*8*q));     %final bitstream to be transmitted.
for i=2:length(str)*8*q
 accbitstream(1,i) = mod(accbitstream(1,i-1) + permbitstream(1,i), 2);
end


%-------UDAIP------------------
accbitstream = permbitstream;
%-------End of UDAIP------------



  bitCount=1;    %number of bits encoded
%------------------------------------------------------------------------------------


[M N] = size(img);
B = 32;                     %B = block size
Mb = floor(M/B);            
Nb = floor(N/B);            %Mb x Nb B sized blocks

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
    newQM = floor((QM * S)/100);
end


newimg=img;     %final encoded image
randBlock = double(zeros(8,8));



for I = 1:Mb
    for J = 1:Nb
        %newimg = img((I-1)*B+1:I*B,(J-1)*B+1:J*B);
        %Randomly selecting an 8x8 block from BxB block
        Sx = randi([1,B-7]);
        Sy = randi([1,B-7]);
        randBlock = double(img((I-1)*B+Sx:(I-1)*B+7+Sx,(J-1)*B+Sy:(J-1)*B+7+Sy));
        %Computing DCT for the randomly selected block
        randBlock = randBlock-128;
        randBlockDCT = dct2(randBlock);
        randBlockQ = (randBlockDCT ./ newQM);
        %randBlockQd = randBlockQ;
        %------------------EMBEDDING-------------------------------------
        zigzag = zeros(1,20);
        zigzag(1)=randBlockQ(1,1);
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
            zigzag(n)=abs(round(randBlockQ(i,j)));
            if(randBlockQ(i,j)>0)
                sgn = '+';
            else
                sgn = '-';
            end
            if(zigzag(n) > 0)                                                   %Checking for threshold
                
                %UDAIP UDAIP-------------------------------
                if(randBlockQ(i,j)>0.5 && randBlockQ(i,j)<0.58)
                    randBlockQ(i,j) = 0.58;
                end
                if(randBlockQ(i,j)<-0.5 && randBlockQ(i,j)>-0.58)
                    randBlockQ(i,j) = -0.58;
                end
                %----End of UDAIP-----%
                if(bitCount<=length(str)*8*q)           
                    if (accbitstream(bitCount)==1)                                  %If bit to be encoded is 1
                        if (mod(zigzag(n),2)==0)                                    %Converting to odd if it is an even reconstruction point
                           % randBlockQ(i,j) = (randBlockDCT(i,j) + newQM(i,j))/newQM(i,j);      %Incrementing to next step size to get odd
                           eval(['randBlockQ(i,j)=(randBlockDCT(i,j)',sgn,'newQM(i,j))/newQM(i,j);']);
                        end
                    elseif  (accbitstream(bitCount)==0)                                 %If bit to be encoded is 0
                        if (mod(zigzag(n),2)~=0)                                        %Converting to even if it is an odd reconstruction point
                            % randBlockQ(i,j) = (randBlockDCT(i,j) + newQM(i,j))/newQM(i,j);      %Incrementing to next step size to get odd                 
                            eval(['randBlockQ(i,j)=(randBlockDCT(i,j)',sgn,'newQM(i,j))/newQM(i,j);']);
                        end
                    end
                bitCount = bitCount+1;
                end
            elseif (zigzag(n)==0)
                randBlockQ(i,j)=zigzag(n);        
            end
        end
        %---------------------------------------------------------------
        
        %Multiplication by JPEG QM and inverse DCT calculation
        irandBlockDCT = randBlockQ .* newQM;
        irandBlock = idct2(irandBlockDCT);
        irandBlock = irandBlock + 128;
        newimg((I-1)*B+Sx:(I-1)*B+7+Sx,(J-1)*B+Sy:(J-1)*B+7+Sy)=uint8(irandBlock);
       
        %Checking loop
%          if(I==1 && J==1)
%             encodedbits = accbitstream(1:bitCount-1);
%             orandBlockQ = randBlockQ;
%             orandBlock = randBlock;
%             save('firstblock', 'orandBlockQ'  , 'bitCount', 'orandBlock', 'irandBlock',  'encodedbits');
%             imwrite(uint8(irandBlock),'test.jpg','JPEG','Mode','lossless');
%          end
        
    end
end
imshow([img newimg]);
imwrite(newimg,'newimage.jpg','Quality',100)%,'JPEG','Mode','lossless');
% clear
% clc
% load('firstblock.mat');
% decode
% isequal(bitstream,encodedbits)
% bc=bitCount;
