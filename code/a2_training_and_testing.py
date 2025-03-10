# Copyright 2020 University of Toronto, all rights reserved

'''Functions related to training and testing.

You don't need anything more than what's been imported here.
'''

import torch
import a2_bleu_score


from tqdm import tqdm


def train_for_epoch(model, dataloader, optimizer, device):
    '''Train an EncoderDecoder for an epoch

    An epoch is one full loop through the training data. This function:

    1. Defines a loss function using :class:`torch.nn.CrossEntropyLoss`,
       keeping track of what id the loss considers "padding"
    2. For every iteration of the `dataloader` (which yields triples
       ``F, F_lens, E``)
       1. Sends ``F`` to the appropriate device via ``F = F.to(device)``. Same
          for ``F_lens`` and ``E``.
       2. Zeros out the model's previous gradient with ``optimizer.zero_grad()``
       3. Calls ``logits = model(F, F_lens, E)`` to determine next-token
          probabilities.
       4. Modifies ``E`` for the loss function, getting rid of a token and
          replacing excess end-of-sequence tokens with padding using
        ``model.get_target_padding_mask()`` and ``torch.masked_fill``
       5. Flattens out the sequence dimension into the batch dimension of both
          ``logits`` and ``E``
       6. Calls ``loss = loss_fn(logits, E)`` to calculate the batch loss
       7. Calls ``loss.backward()`` to backpropagate gradients through
          ``model``
       8. Calls ``optimizer.step()`` to update model parameters
    3. Returns the average loss over sequences

    Parameters
    ----------
    model : EncoderDecoder
        The model we're training.
    dataloader : HansardDataLoader
        Serves up batches of data.
    device : torch.device
        A torch device, like 'cpu' or 'cuda'. Where to perform computations.
    optimizer : torch.optim.Optimizer
        Implements some algorithm for updating parameters using gradient
        calculations.

    Returns
    -------
    avg_loss : float
        The total loss divided by the total numer of sequence
    '''
    # If you want, instead of looping through your dataloader as
    # for ... in dataloader: ...
    # you can wrap dataloader with "tqdm":
    # for ... in tqdm(dataloader): ...
    # This will update a progress bar on every iteration that it prints
    # to stdout. It's a good gauge for how long the rest of the epoch
    # will take. This is entirely optional - we won't grade you differently
    # either way.
    # If you are running into CUDA memory errors part way through training,
    # try "del F, F_lens, E, logits, loss" at the end of each iteration of
    # the loop.
    criterion = torch.nn.CrossEntropyLoss(ignore_index = model.target_vocab_size)
    total_loss = 0
    iteration = 0
    criterion = criterion.to(device)
    model = model.to(device)
    for F, F_lens, E in tqdm(dataloader):
        #Sends ``F`` to the appropriate device via ``F = F.to(device)``. Same
        #for ``F_lens`` and ``E``.
        F = F.to(device)
        F_lens = F_lens.to(device)
        E = E.to(device)

        #Zeros out the model's previous gradient with ``optimizer.zero_grad()``
        optimizer.zero_grad()
        #Calls ``logits = model(F, F_lens, E)`` to determine next-token probabilities.
        logits = model(F, F_lens, E)
        #logits = logits.to(device)

        #Modifies ``E`` for the loss function, getting rid of a token and
        #replacing excess end-of-sequence tokens with padding using
        #``model.get_target_padding_mask()`` and ``torch.masked_fill
        padding_mask = model.get_target_padding_mask(E)
        E = E.masked_fill(padding_mask, model.target_vocab_size)

        #Flattens out the sequence dimension into the batch dimension of both
        #``logits`` and ``E``
        logits = torch.flatten(logits, 0, -2)
        E = torch.flatten(E[1:], 0)
        #E = torch.flatten(E, 0)

        loss = criterion(logits, E)
        iteration += 1
        total_loss += loss
        loss.backward()
        optimizer.step()
        del F, F_lens, E, logits, loss
    avg_loss = total_loss / iteration
    #print("Average Loss: " + str(avg_loss))
    #print("Number of batches: " + str(iteration))
    return avg_loss







def compute_batch_total_bleu(E_ref, E_cand, target_sos, target_eos):
    '''Compute the total BLEU score over elements in a batch

    Parameters
    ----------
    E_ref : torch.LongTensor
        A batch of reference transcripts of shape ``(T, N)``, including
        start-of-sequence tags and right-padded with end-of-sequence tags.
    E_cand : torch.LongTensor
        A batch of candidate transcripts of shape ``(T', N)``, also including
        start-of-sequence and end-of-sequence tags.
    target_sos : int
        The ID of the start-of-sequence tag in the target vocabulary.
    target_eos : int
        The ID of the end-of-sequence tag in the target vocabulary.

    Returns
    -------
    total_bleu : float
        The sum total BLEU score for across all elements in the batch. Use
        n-gram precision 4.
    '''
    # you can use E_ref.tolist() to convert the LongTensor to a python list
    # of numbers
    total_bleu = 0
    reference = E_ref.tolist()
    candidate = E_cand.tolist()
    for i, j in zip(range(len(reference)), range(len(candidate))):
        ref = reference[i]
        cand = candidate[j]
        str_ref = []
        str_cand = []
        str_target_eos = str(target_eos)
        str_target_sos = str(target_sos)
        for s in ref:
            str_ref.append(str(s))
        for s in cand:
            str_cand.append(str(s))
        #print("cand Length" + str(len(str_cand)))
        #print(str_cand)
        #print("ref Length" + str(len(str_ref)))
        if str_target_sos in str_ref:
            str_ref.remove(str_target_sos)

        if str_target_eos in str_ref:
            str_ref.remove(str_target_eos)

        if str_target_sos in str_cand:
            str_cand.remove(str_target_sos)

        if str_target_eos in str_cand:
            str_cand.remove(str_target_eos)
        bleu_score = a2_bleu_score.BLEU_score(str_ref, str_cand, 4)
        total_bleu += bleu_score
    return total_bleu



def compute_average_bleu_over_dataset(
        model, dataloader, target_sos, target_eos, device):
    '''Determine the average BLEU score across sequences

    This function computes the average BLEU score across all sequences in
    a single loop through the `dataloader`.

    1. For every iteration of the `dataloader` (which yields triples
       ``F, F_lens, E_ref``):
       1. Sends ``F`` to the appropriate device via ``F = F.to(device)``. Same
          for ``F_lens``. No need for ``E_cand``, since it will always be
          compared on the CPU.
       2. Performs a beam search by calling ``b_1 = model(F, F_lens)``
       3. Extracts the top path per beam as ``E_cand = b_1[..., 0]``
       4. Computes the total BLEU score of the batch using
          :func:`compute_batch_total_bleu`
    2. Returns the average per-sequence BLEU score

    Parameters
    ----------
    model : EncoderDecoder
        The model we're testing.
    dataloader : HansardDataLoader
        Serves up batches of data.
    target_sos : int
        The ID of the start-of-sequence tag in the target vocabulary.
    target_eos : int
        The ID of the end-of-sequence tag in the target vocabulary.

    Returns
    -------
    avg_bleu : float
        The total BLEU score summed over all sequences divided by the number of
        sequences
    '''
    model = model.to(device)
    total_bleu_score = 0
    iterations = 0
    for F, F_lens, E_ref in tqdm(dataloader):
        F = F.to(device)
        F_lens = F_lens.to(device)
        b_1 = model(F, F_lens)
        E_cand = b_1[..., 0]
        total_bleu_score += compute_batch_total_bleu(E_ref, E_cand, target_sos, target_eos)
        iterations += 1
    avg_bleu = total_bleu_score / (iterations * 100)
    #print("Average Bleu Score: " + str(avg_bleu))
    #print("Number of batches: " + str(iterations))
    return avg_bleu


