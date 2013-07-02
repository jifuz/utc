I=prior;

figure
subplot(1,2,1); 
imshow(I); title('Original Prior');

H = fspecial('motion',30,0);
H = fspecial('gaussian', [1, 50], 100);

Itemp=I;
Itemp(Itemp==0)=-1;
Itemp(Itemp>0)=0;
Itemp=-Itemp;
MotionBlur = imfilter(Itemp,H,'replicate');
MotionBlur=max(MotionBlur(:))-MotionBlur;
MotionBlur(Itemp>0)=0;
subplot(1,2,2); 
imshow(MotionBlur);title('Blurred Prior');



original = imread('cameraman.tif');

original=zeros(300,300);
original(100:200,100:200)=1;
temp=original(100:200,100:200);
PSF = fspecial('gaussian',40,10);
edgesTapered = edgetaper(temp,PSF);
original(100:200,100:200)=edgesTapered;
figure, imshow(original,[]);
figure, imshow(edgesTapered,[]);