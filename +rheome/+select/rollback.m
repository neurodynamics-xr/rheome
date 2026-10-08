function n = rollback(db, txnId)
% SELECT.ROLLBACK  Remove a transaction's rows and its txn row, atomically.
%
%   n = rheome.select.rollback(db, txnId)     -> rows removed
%
% Labels and measurements share one transaction log, so one call undoes whichever the
% transaction wrote; a write is never half-undone.
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        txnId (1,1) double
    end
    [L, T, M] = sel_labelio(db);
    if ~any(T.txn_id == txnId)
        error('select:rollback:txn', 'No transaction %d in %s.', txnId, db.labelFile);
    end
    keepL = L.txn_id ~= txnId;  keepM = M.txn_id ~= txnId;
    n = nnz(~keepL) + nnz(~keepM);
    sel_labelio(db, L(keepL, :), T(T.txn_id ~= txnId, :), M(keepM, :));
end
% Author: Diellor Basha, 2026
