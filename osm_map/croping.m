original = imread('itu_sat.jpg');

org_size = size(original);
width  = org_size(2);
heigth = org_size(1);

rect_w = 200;
rect_h = 200;

w_slid = 200;
h_slid = 100;

w_iter = fix((width - rect_w) / w_slid);
h_iter = fix((heigth- rect_h) / h_slid);

counter = 0;
for i=1:w_iter
    for j=1:h_iter 

        [J,rect] = imcrop(original,[(i-1)*w_slid (j-1)*h_slid rect_w rect_h]);
        imwrite(J,['crop_imgs\croped_', num2str(counter), '.jpg']);

        counter = counter + 1;
    end
end

