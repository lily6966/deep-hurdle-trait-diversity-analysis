import torch
import torch.nn.functional as F 
import torch.nn as nn


class VAE(nn.Module):
    def __init__(self, args):
        super(VAE, self).__init__()
        self.args = args
        # feature layers
        input_dim = args.input_dim #+ args.meta_offset
        self.fx1 = nn.Linear(input_dim, 256)
        self.fx2 = nn.Linear(256, 512)
        self.fx3 = nn.Linear(512, 256)
        self.fx_mu = nn.Linear(256, args.latent_dim)
        self.fx_logvar = nn.Linear(256, args.latent_dim)

        self.emb_size = args.emb_size

        self.fd_x1 = nn.Linear(input_dim + args.latent_dim, 512)
        self.fd_x2 = torch.nn.Sequential(
            nn.Linear(512, self.emb_size)
        )
        self.feat_mp_mu = nn.Linear(self.emb_size, args.label_dim)

        self.recon = torch.nn.Sequential(
            nn.Linear(args.latent_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, input_dim)
        )

        self.label_recon = torch.nn.Sequential(
            nn.Linear(args.latent_dim, 512),
            nn.ReLU(),
            nn.Linear(512, self.emb_size),
            nn.LeakyReLU()
        )

        # label layers
        self.fe0 = nn.Linear(args.label_dim, self.emb_size)
        self.fe1 = nn.Linear(self.emb_size, 512)
        self.fe2 = nn.Linear(512, 256)
        self.fe_mu = nn.Linear(256, args.latent_dim)
        self.fe_logvar = nn.Linear(256, args.latent_dim)

        self.fd1 = self.fd_x1
        self.fd2 = self.fd_x2
        self.label_mp_mu = self.feat_mp_mu

        self.bias = nn.Parameter(torch.zeros(args.label_dim))

        assert id(self.fd_x1) == id(self.fd1)
        assert id(self.fd_x2) == id(self.fd2)

        self.dropout = nn.Dropout(p=args.keep_prob)
        self.scale_coeff = args.scale_coeff

    def label_encode(self, x):
        h0 = self.dropout(F.relu(self.fe0(x)))  # [label_dim, emb_size]
        h1 = self.dropout(F.relu(self.fe1(h0)))  # [label_dim, 512]
        h2 = self.dropout(F.relu(self.fe2(h1)))  # [label_dim, 256]
        mu = self.fe_mu(h2) * self.scale_coeff  # [label_dim, latent_dim]
        logvar = self.fe_logvar(h2) * self.scale_coeff  # [label_dim, latent_dim]

        fe_output = {
            'fe_mu': mu,
            'fe_logvar': logvar
        }
        return fe_output

    def feat_encode(self, x):
        h1 = self.dropout(F.relu(self.fx1(x)))
        h2 = self.dropout(F.relu(self.fx2(h1)))
        h3 = self.dropout(F.relu(self.fx3(h2)))
        mu = self.fx_mu(h3) * self.scale_coeff  # [bs, latent_dim]
        logvar = self.fx_logvar(h3) * self.scale_coeff
        fx_output = {
            'fx_mu': mu,
            'fx_logvar': logvar
        }
        return fx_output

    def label_reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def feat_reparameterize(self, mu, logvar, coeff=1.0):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def label_decode(self, z):
        d1 = F.relu(self.fd1(z))
        d2 = F.leaky_relu(self.fd2(d1))
        d3 = F.normalize(d2, dim=1)
        return d3

    def feat_decode(self, z):
        d1 = F.relu(self.fd_x1(z))
        d2 = F.leaky_relu(self.fd_x2(d1))
        d3 = F.normalize(d2, dim=1)
        return d3

    def label_forward(self, x, feat):  # x is label
        if self.args.reg == "gmvae":
            n_label = x.shape[1]  # label_dim
            all_labels = torch.eye(n_label).to(x.device)  # [label_dim, label_dim]
            fe_output = self.label_encode(all_labels)  # map each label to a Gaussian mixture.
            # print(fe_output['fe_mu'].shape) # [label_dim, latent_dim]
        else:
            fe_output = self.label_encode(x)
        mu = fe_output['fe_mu']
        logvar = fe_output['fe_logvar']

        if self.args.reg == "wae" or not self.training:
            if self.args.reg == "gmvae":
                z = torch.matmul(x, mu) / (x.sum(1, keepdim=True)+1e-8)
            else:
                z = mu
        else:
            if self.args.reg == "gmvae":
                z = torch.matmul(x, mu) / (x.sum(1, keepdim=True)+1e-8)  # mu of Gaussian Mixture
                # t = x.sum(1, keepdim=True) # [bs, 1]
                # print(x.shape, mu.shape) [bs, label_dim], [label_dim, latent_dim]
                # print(z.shape) # [bs, latent_dim]
            else:
                z = self.label_reparameterize(mu, logvar)
        label_emb = self.label_decode(torch.cat((feat, z), 1))
        # print(label_emb.shape) # [bs, emb_size]
        single_label_emb = F.normalize(self.label_recon(mu), dim=1)  # [label_dim, emb_size]. Seems not useful

        fe_output['label_emb'] = label_emb
        fe_output['single_label_emb'] = single_label_emb
        return fe_output

    def feat_forward(self, x):
        fx_output = self.feat_encode(x)
        mu = fx_output['fx_mu']  # [bs, latent_dim]
        logvar = fx_output['fx_logvar']  # [bs, latent_dim]

        if self.args.reg == "wae" or not self.training:
            if self.args.test_sample:
                z = self.feat_reparameterize(mu, logvar)
                z2 = self.feat_reparameterize(mu, logvar)
            else:
                z = mu
                z2 = mu
        else:
            z = self.feat_reparameterize(mu, logvar)  # [bs, latent_dim]
            z2 = self.feat_reparameterize(mu, logvar)  # [bs, latent_dim]

        feat_emb = self.feat_decode(torch.cat((x, z), 1))  # [bs, emb_size]
        feat_emb2 = self.feat_decode(torch.cat((x, z2), 1))  # [bs, emb_size]
        fx_output['feat_emb'] = feat_emb
        fx_output['feat_emb2'] = feat_emb2

        feat_recon = self.recon(z)
        fx_output['feat_recon'] = feat_recon
        return fx_output

    def forward(self, input_feat, input_label_count, input_label_binary):
        # input_feat, input_label_count, input_label_binary
        '''
        if self.args.mode == 'classification':
            fe_output = self.label_forward(input_label_binary, input_feat)
        elif self.args.mode == 'regression':
            fe_output = self.label_forward(input_label_count, input_feat)
        else:
            print("mode error!!!")
            exit()
        '''
        fe_output = self.label_forward(input_label_binary, input_feat)
        label_emb, single_label_emb = fe_output['label_emb'], fe_output[
            'single_label_emb']  # [bs, emb_size], [label_dim, emb_size]
        fx_output = self.feat_forward(input_feat)
        feat_emb, feat_emb2 = fx_output['feat_emb'], fx_output['feat_emb2']  # [bs, emb_size], [bs, emb_size]
        embs = self.fe0.weight
        # print(embs.shape) # [emb_size, label_dim]

        label_out = torch.matmul(label_emb, embs)  # [bs, emb_size] * [emb_size, label_dim] = [bs, label_dim]
        single_label_out = torch.matmul(single_label_emb, embs)  # [label_dim, label_dim]

        feat_out = torch.matmul(feat_emb, embs)  # [bs, label_dim]
        feat_out2 = torch.matmul(feat_emb2, embs)  # [bs, label_dim]

        fe_output.update(fx_output)
        output = fe_output
        output['embs'] = embs
        output['label_out'] = label_out
        output['single_label_out'] = single_label_out
        output['feat_out'] = feat_out
        output['feat_out2'] = feat_out2
        output['feat'] = input_feat

        return output


def compute_loss(input_label_binary, input_label_count, output, criterion, args=None, epoch=0,
                 class_weights=None, mode='classification', pred_binary=None):

    def log_sum_exp(x, mask):
        max_x = torch.max(x, 1)[0]
        new_x = x - max_x.unsqueeze(1).expand_as(x)
        return max_x + (new_x.exp().sum(1)).log()

    def log_mean_exp(x, mask):
        return log_sum_exp(x, mask) - torch.log(mask.sum(1)+1e-8)

    def log_normal(x, m, v):
        log_prob = (-0.5 * (torch.log(v+1e-8) + (x - m).pow(2) / (v+1e-8))).sum(-1)
        return log_prob

    def log_normal_mixture(z, m, v, mask=None):
        m = m.unsqueeze(0).expand(z.shape[0], -1, -1)
        v = v.unsqueeze(0).expand(z.shape[0], -1, -1)
        batch, mix, dim = m.size()
        z = z.view(batch, 1, dim).expand(batch, mix, dim)
        indiv_log_prob = log_normal(z, m, v) + torch.ones_like(mask) * (-1e6) * (1. - mask)
        log_prob = log_mean_exp(indiv_log_prob, mask)
        return log_prob

    fe_out, fe_mu, fe_logvar, label_emb = output['label_out'], output['fe_mu'], output['fe_logvar'], output[
        'label_emb']
    fx_out, fx_mu, fx_logvar, feat_emb = output['feat_out'], output['fx_mu'], output['fx_logvar'], output[
        'feat_emb']
    fx_out2, single_label_out = output['feat_out2'], output['single_label_out']
    embs = output['embs']

    feat_recon_loss = 0.

    fe_sample = torch.matmul(input_label_binary, fe_mu) / (input_label_binary.sum(1, keepdim=True)+1e-8)

    std = torch.exp(0.5 * fx_logvar)
    eps = torch.randn_like(std)
    fx_sample = fx_mu + eps * std
    fx_var = torch.exp(fx_logvar)
    fe_var = torch.exp(fe_logvar)
    kl_loss = (log_normal(fx_sample, fx_mu, fx_var) - log_normal_mixture(fx_sample, fe_mu, fe_var,
                                                                             input_label_binary)).mean()

    def supconloss(label_emb, feat_emb, embs, temp=1.0, sample_wise=False):
        if sample_wise:
            loss_func = SupConLoss(temperature=0.1)
            return loss_func(torch.stack([label_emb, feat_emb], dim=1), input_label_binary.float())

        features = torch.cat((label_emb, feat_emb))
        labels = torch.cat((input_label_binary, input_label_binary)).float()
        n_label = labels.shape[1]
        emb_labels = torch.eye(n_label).to(feat_emb.device)
        mask = torch.matmul(labels, emb_labels)

        anchor_dot_contrast = torch.div(torch.matmul(features, embs),temp+1e-8)
        logits_max, _ = torch.max(anchor_dot_contrast, dim=1, keepdim=True)
        logits = anchor_dot_contrast - logits_max.detach()

        exp_logits = torch.exp(logits)
        log_prob = logits - torch.log(exp_logits.sum(1, keepdim=True)+1e-8)

        mean_log_prob_pos = (mask * log_prob).sum(1) / (mask.sum(1)+1e-8)
        mean_log_prob_neg = ((1.0 - mask) * log_prob).sum(1) / ((1.0 - mask).sum(1)+1e-8)

        loss = - mean_log_prob_pos
        loss = loss.mean()

        return loss

    def compute_BCE(E):
        # compute negative log likelihood (BCE loss) for each sample point
        sample_nll = -(torch.log(E) * input_label_binary + torch.log(1 - E) * (1 - input_label_binary))
        logprob = -torch.sum(sample_nll, dim=2)

        # the following computation is designed to avoid the float overflow (log_sum_exp trick)
        maxlogprob = torch.max(logprob, dim=0)[0]
        Eprob = torch.mean(torch.exp(logprob - maxlogprob), axis=0)
        nll_loss = torch.mean(-torch.log(Eprob) - maxlogprob)
        # c_loss = build_multi_classification_loss(E, input_label)
        return nll_loss

    if mode == 'classification':
        pred_e = torch.sigmoid(fe_out)
        pred_x = torch.sigmoid(fx_out)
        pred_x2 = torch.sigmoid(fx_out2)
        #nll_loss = criterion(pred_e, input_label_binary)
        #nll_loss_x = criterion(pred_x, input_label_binary)
        #nll_loss_x2 = criterion(pred_x2, input_label_binary)

        nll_loss = compute_BCE(pred_e.unsqueeze(0))
        nll_loss_x = compute_BCE(pred_x.unsqueeze(0))
        nll_loss_x2 = compute_BCE(pred_x2.unsqueeze(0))

    elif mode == 'regression':
        l1_loss = torch.nn.L1Loss(reduction='none')
        pred_e = torch.exp(fe_out)
        pred_x = torch.exp(fx_out)
        pred_x2 = torch.exp(fx_out2)
        nll_loss = torch.mean(criterion(pred_e,input_label_count)*input_label_binary) + 0.01*torch.mean(l1_loss(pred_e,input_label_count))
        nll_loss_x = torch.mean(criterion(pred_x,input_label_count)*input_label_binary) + 0.01*torch.mean(l1_loss(pred_x,input_label_count))
        nll_loss_x2 = torch.mean(criterion(pred_x2,input_label_count)*input_label_binary) + 0.01*torch.mean(l1_loss(pred_x2,input_label_count))

    cpc_loss = supconloss(label_emb, feat_emb, embs, sample_wise=False)
    total_loss = (nll_loss + nll_loss_x + nll_loss_x2) * 10 + kl_loss * 6. + 1*(cpc_loss)  # + latent_cpc_loss

    return total_loss, nll_loss, nll_loss_x, kl_loss, cpc_loss, pred_e, pred_x
