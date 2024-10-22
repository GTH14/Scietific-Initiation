%% FUNÇÕES
function p_total = Calc_pressao_estacionario(x, y, z, transdutor)
    x_t = transdutor.posx;
    y_t = transdutor.posy;
    z_t = transdutor.posz;
    A = transdutor.A;
    normal_v = transdutor.normal;
%     normal_v2 = transdutor.normal2;
    omega = transdutor.omega;
    a = transdutor.a;
    k = transdutor.k;
    phase_v = transdutor.phase;
    n1 = length(x_t);
    s = size(x);
    r = zeros([s,n1])+0j;
    theta = zeros([s,n1])+0j;
    Pax = zeros([s,n1])+0j;
    v = zeros([s,n1])+0j;
    H = zeros([s,n1])+0j;
    p = zeros([s,n1])+0j;
    t1 = 0;
    for i = 1:n1
        
        r(:,:,i) = calc_r(x_t(i), y_t(i), z_t(i), x, y,z);
%         if i<=n1/2
            theta(:,:,i) = calc_th(x_t(i), y_t(i), z_t(i), x, y, z, normal_v, r(:,:,i));
%         else
%             theta(:,:,i) = calc_th(x_t(i), y_t(i), z_t(i), x, y, z, normal_v2, r(:,:,i));
%         end
        Pax(:,:,i) = calc_Pax(A, r(:,:,i));
        v(:,:,i) = calc_v(k, a, theta(:,:,i));
        H(:,:,i) = calc_H(v(:,:,i));
        p(:,:,i) = calc_p(Pax(:,:,i), H(:,:,i), omega, k, t1, r(:,:,i), phase_v(i));
    end
    p_total = sum(p,3);
end
