function W = plotMorletWavelet(nx,x0,f)

sigma = 1/f;
X = 1:nx
W = exp(-(X-x0).^2/sigma^2).*exp(2i*pi*f*(X-x0));

figure();
plot(X,abs(W),'k',X,real(W),'b',X,imag(W),'r');
end