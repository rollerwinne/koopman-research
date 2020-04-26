function z=rate(Z,r)
zmax=max(Z(:));
zmin=min(Z(:));
z=zmin+r*(zmax-zmin);
end