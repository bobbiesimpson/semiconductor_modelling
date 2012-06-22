function L2pnJunctionPlots( DOF, L2psi, L2p, L2n)

figure(1)
grid on
hold on
loglog(DOF, L2psi, 'ko--');
hold off
xlabel('DOF')
ylabel('L2 potential norm')


figure(2)
grid on
hold on
loglog(DOF, L2p, 'ko--');
hold off
xlabel('DOF')
ylabel('L2 hole concentration norm')


figure(3)
grid on
hold on
loglog(DOF, L2n, 'ko--');
hold off
xlabel('DOF')
ylabel('L2 electron concentration norm')

psiExport = [DOF' L2psi'];
pExport = [DOF' L2p'];
nExport = [DOF' L2n'];

save 'datFiles/L2psi.dat' 'psiExport' '-ASCII'
save 'datFiles/L2p.dat' 'pExport' '-ASCII'
save 'datFiles/L2n.dat' 'nExport' '-ASCII'

end

