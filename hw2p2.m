% HW2 p2
 im1 = iradon(sinogramD', [0:359]);
 im2 = iradon(sinogramD(1:2:end,:)', [0:2:359]);
 im3 = iradon(sinogramD(1:4:end,:)', [0:4:359]);
 im4 = iradon(sinogramD(1:6:end,:)', [0:6:359]);
 
 figure; imshow(im1, []); colorbar;
 export_fig('hw2p2a1','-pdf');
 
 figure; imshow(im2, []); colorbar;
 export_fig('hw2p2a2','-pdf');
 figure; imshow(im3, []); colorbar;
 export_fig('hw2p2a3','-pdf');
 figure; imshow(im4, []); colorbar;
 export_fig('hw2p2a4','-pdf');
 
