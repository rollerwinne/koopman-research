% function [KV,KD]=Koopman_svd_draw(K,L,r)
% [U,E,V]=svd(K);
% U=U(:,1:r);
% E=E(1:r,1:r);
% V=V(:,1:r);
% K_new=E*V';
% L_new=U'*L;
% A=L_new*pinv(K_new);
% [KV,KD]=eig(A);
% KV=U*KV;
% end

function [KV,KD]=Koopman_svd_draw(K,L,r)
[U,E,V]=svd(K);
U=U(:,1:r);
E=E(1:r,1:r);
V=V(:,1:r);
K_new=E*V';
L_new=U'*L;
A=pinv(K_new)*L_new;
[KV,KD]=eig(A);
KV=K_new*KV;
KV=U*KV;
end