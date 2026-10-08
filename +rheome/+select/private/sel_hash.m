function h = sel_hash(s)
% SEL_HASH  MD5 of a string, lower-case hex: the content hash a transaction is idempotent on.
%
% One definition for labels and measurements, because two transactions that wrote the same
% rows must collide, and a hash that differed by table would let the same write land twice.
%
% Author: Diellor Basha, 2026
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(s));
    d = typecast(md.digest(), 'uint8');
    h = lower(reshape(dec2hex(d, 2)', 1, []));
end
% Author: Diellor Basha, 2026
