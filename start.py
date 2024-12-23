import numpy

frame_bit_stream: str, curr_frame_resd: np.ndarray = RPE_frame_coder(
        s0: np.ndarray,
        prev_frame_resd: np.ndarray
        )

# s0: ένα numpy array, που αναπαριστά ένα διάνυσμα από 160 δείγματα φωνής
# prev_frame_resd: οι τιμές της ακολουθίας d′(n) μήκους 160 δειγμάτων, του προηγούμενου frame
# frame_bit_stream: τα 260 bits που αναπαριστούν το τρέχον frame.
#   Σημείωση: η αναπαράσταση του bitstream μπορεί για ευκολία να είναι ένα απλό string από τους χαρακτήρες
#   '0','1', είτε να χρησιμοποιεί τη βιβλιοθήκη bitstream1, είτε κάποια άλλη υλοποίηση. Παρόλα αυτά, κάθε
#   επιλογή πρέπει να περιγράφεται επαρκώς στην αναφορά σας.
# curr_frame_resd: οι τιμές της ακολουθίας d′(n) μήκους 160 δειγμάτων, του τρέχοντος frame. 
#
#
# Χρειάζεται να φτιάξουμε ένα wrapper πρόγραμμα που διαβάζει/γράφει ανά 160 δείγματα αρχεία τύπου .wav
#
#
# # # # # #
# LEVEL 1 #
# # # # # #
#
# LARc: np.ndarray, curr_frame_st_resd: np.ndarray = RPE_frame_st_coder(
#       s0: np.ndarray
#       )
#
# s0: np.ndarray = RPE_frame_st_decoder(
#       LARc: np.ndarray,
#       curr_frame_st_resd: np.ndarray
#       )


