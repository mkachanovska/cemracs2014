function [err1, err2, err2l2, epw1, epw2]=compute_error(Kstiff, Mmass, e1, e2, e1c, e2c)
err1p=e1-e1c;
err1=sqrt(err1p'*Mmass*err1p);
err2p=e2-e2c;
err2=sqrt(err2p'*Mmass*err2p+err2p'*Kstiff*err2p);

err2l2=sqrt(err2p'*Mmass*err2p);

epw1=norm(err1p);
epw2=norm(err2p);
end
