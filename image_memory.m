
img = randi([-1,1],10,10);
img(img==0) = 1;

imwrite(img,'pic.jpg');

imshow('pic.jpg');
