function ueu = export_phase(phases)
    fid = fopen( 'Teste_1.txt', 'wt' );
    n = length(phase_v);
    for i = 1:n
        formatSpec = 'phase'+ '%4.2f meters or %8.3f mm\n';
        fprint('Teste_1.txt',
    end

end
