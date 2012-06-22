function  plotCurrentVariables( meshcoords, currentPsi, current_n, current_p )
% plot the current variables

figure(1);
hold on
plot(meshcoords, currentPsi,'ks--')
grid on;
xlabel('x (cm)');
ylabel('potential (V)');
hold off

export2 = [meshcoords' currentPsi];
save 'dat_files/psi.dat' export2 -ASCII;

figure(2);
hold on
plot(meshcoords, current_n, 'ks--');
xlabel('x (cm)');
ylabel('electron concentration (cm^{-3})');
grid on;
hold off

export2 = [meshcoords' current_n];
save 'dat_files/n.dat' export2 -ASCII;

figure(3);
hold on
plot(meshcoords, current_p, 'ks--');
xlabel('x (cm)');
ylabel('hole concentration (cm^{-3})');
grid on;
hold off

export2 = [meshcoords' current_p];
save 'dat_files/p.dat' export2 -ASCII;

end

