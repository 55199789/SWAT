/**
*
* MIT License
*
* Copyright (c) Open Enclave SDK contributors.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE
*
*/

#include "common/openssl_utility.h"
#include <errno.h>

sgx_status_t generate_certificate_and_pkey(X509*& certificate, EVP_PKEY*& pkey)
{
    quote3_error_t qresult = SGX_QL_SUCCESS;
    sgx_status_t result = SGX_ERROR_UNEXPECTED;
    uint8_t* output_certificate = NULL;
    size_t output_certificate_size = 0;
    uint8_t* private_key_buffer = nullptr;
    size_t private_key_buffer_size = 0;
    uint8_t* public_key_buffer = nullptr;
    size_t public_key_buffer_size = 0;
    const uint8_t* certificate_buffer_ptr = nullptr;
    BIO* mem = nullptr;
    int key_type = RSA_TYPE;

    if (key_type) {
        t_print(" generating keys by EC P-384\n");
    }
    else
    {
        t_print(" generating keys by RSA 3072\n");
    }
    result = generate_key_pair(
        key_type, &public_key_buffer,
        &public_key_buffer_size,
        &private_key_buffer,
        &private_key_buffer_size);
    if (result != SGX_SUCCESS)
    {
        t_print(" failed to generate RSA key pair\n");
        goto done;
    }

    // t_print("certificate addr = %p\n", certificate_subject_name);
    // t_print("public_key_buf_size:[%ld], addr = %p\n", public_key_buffer_size, public_key_buffer);
    // t_print("%s\n", public_key_buffer);
    // t_print("private_key_buf_size:[%ld], addr = %p\n", private_key_buffer_size, private_key_buffer);
    // t_print("%s\n", private_key_buffer);

    qresult = tee_get_certificate_with_evidence(
        certificate_subject_name,
        private_key_buffer,
        private_key_buffer_size,
        public_key_buffer,
        public_key_buffer_size,
        &output_certificate,
        &output_certificate_size);

    if (qresult != SGX_QL_SUCCESS || output_certificate == nullptr)
    {
        if (output_certificate == nullptr)
            t_print(" null certificate\n");
        p_sgx_tls_qe_err_msg(qresult);
        goto done;
    }

    // temporary buffer required as if d2i_x509 call is successful
    // certificate_buffer_ptr is incremented to the byte following the parsed
    // data. sending certificate_buffer_ptr as argument will keep
    // output_certificate pointer undisturbed.

    certificate_buffer_ptr = output_certificate;

    if ((certificate = d2i_X509(
        nullptr,
        &certificate_buffer_ptr,
        (long)output_certificate_size)) == nullptr)
    {
        t_print("Failed to convert DER format certificate to X509 structure\n");
        goto done;
    }
    mem = BIO_new_mem_buf((void*)private_key_buffer, -1);
    if (!mem)
    {
        t_print("Failed to convert private key buf into BIO_mem\n");
        goto done;
    }
    if ((pkey = PEM_read_bio_PrivateKey(mem, nullptr, 0, nullptr)) == nullptr)
    {
        t_print("Failed to convert private key buffer into EVP_KEY format\n");
        goto done;
    }

    result = SGX_SUCCESS;
done:
    if (private_key_buffer)
        free(private_key_buffer);
    if (public_key_buffer)
        free(public_key_buffer);
    certificate_buffer_ptr = nullptr;

    if (mem)
        BIO_free(mem);
    if (output_certificate)
        tee_free_certificate(output_certificate);
    return result;
}

sgx_status_t load_tls_certificates_and_keys(
    SSL_CTX* ctx,
    X509*& certificate,
    EVP_PKEY*& pkey)
{
    sgx_status_t result = SGX_ERROR_UNEXPECTED;

    if (generate_certificate_and_pkey(certificate, pkey) != SGX_SUCCESS)
    {
        t_print("Cannot generate certificate and pkey\n");
        goto exit;
    }

    if (certificate == nullptr)
    {
        t_print("null cert\n");
        goto exit;
    }

    if (!SSL_CTX_use_certificate(ctx, certificate))
    {
        t_print("Cannot load certificate on the server\n");
        goto exit;
    }

    if (!SSL_CTX_use_PrivateKey(ctx, pkey))
    {
        t_print("Cannot load private key on the server\n");
        goto exit;
    }

    /* verify private key */
    if (!SSL_CTX_check_private_key(ctx))
    {
        t_print("Private key does not match the public certificate\n");
        goto exit;
    }
    result = SGX_SUCCESS;
exit:
    return result;
}

sgx_status_t initalize_ssl_context(SSL_CONF_CTX*& ssl_conf_ctx, SSL_CTX*& ctx)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    // Configure the SSL context based on Open Enclave's security guidance.
    const char* cipher_list_tlsv12_below =
        "ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-ECDSA-AES256-GCM-SHA384:ECDHE-RSA-"
        "AES128-GCM-SHA256:"
        "ECDHE-RSA-AES256-GCM-SHA384:ECDHE-ECDSA-AES128-SHA256:ECDHE-ECDSA-"
        "AES256-SHA384:"
        "ECDHE-RSA-AES128-SHA256:ECDHE-RSA-AES256-SHA384";
    const char* cipher_list_tlsv13 =
        "TLS13-AES-256-GCM-SHA384:TLS13-AES-128-GCM-SHA256";
    const char* supported_curves = "P-521:P-384:P-256";

    SSL_CONF_CTX_set_ssl_ctx(ssl_conf_ctx, ctx);
    SSL_CONF_CTX_set_flags(
        ssl_conf_ctx,
        SSL_CONF_FLAG_FILE | SSL_CONF_FLAG_SERVER | SSL_CONF_FLAG_CLIENT);
    int ssl_conf_return_value = -1;
    if ((ssl_conf_return_value =
        SSL_CONF_cmd(ssl_conf_ctx, "MinProtocol", "TLSv1.2")) < 0)
    {
        t_print(
            "Setting MinProtocol for ssl context configuration failed with "
            "error %d \n",
            ssl_conf_return_value);
        goto exit;
    }
    if ((ssl_conf_return_value =
        SSL_CONF_cmd(ssl_conf_ctx, "MaxProtocol", "TLSv1.3")) < 0)
    {
        t_print(
            "Setting MaxProtocol for ssl context configuration failed with "
            "error %d \n",
            ssl_conf_return_value);
        goto exit;
    }
    if ((ssl_conf_return_value = SSL_CONF_cmd(
        ssl_conf_ctx, "CipherString", cipher_list_tlsv12_below)) < 0)
    {
        t_print(
            "Setting CipherString for ssl context configuration failed with "
            "error %d \n",
            ssl_conf_return_value);
        goto exit;
    }
    if ((ssl_conf_return_value = SSL_CONF_cmd(
        ssl_conf_ctx, "Ciphersuites", cipher_list_tlsv13)) < 0)
    {
        t_print(
            "Setting Ciphersuites for ssl context configuration failed with "
            "error %d \n",
            ssl_conf_return_value);
        goto exit;
    }
    if ((ssl_conf_return_value =
        SSL_CONF_cmd(ssl_conf_ctx, "Curves", supported_curves)) < 0)
    {
        t_print(
            "Setting Curves for ssl context configuration failed with error %d "
            "\n",
            ssl_conf_return_value);
        goto exit;
    }
    if (!SSL_CONF_CTX_finish(ssl_conf_ctx))
    {
        t_print("Error finishing ssl context configuration \n");
        goto exit;
    }
    ret = SGX_SUCCESS;
exit:
    return ret;
}

int read_from_session_peer(
    SSL*& ssl_session,
    void* buffer,
    size_t bufferSize)
{
    int ret = 0;
    int bytesRead = 0;
    size_t accBytesRead = 0;
    do
    {
        bytesRead = SSL_read(ssl_session, buffer, bufferSize);

        if (bytesRead <= 0)
        {
            int error = SSL_get_error(ssl_session, bytesRead);
            if (error == SSL_ERROR_WANT_READ)
                continue;

            t_print("Failed! SSL_read returned error=%d, errno = %s\n", error, strerror(errno));
            return accBytesRead;
        }

        accBytesRead += bytesRead;
        buffer += bytesRead;
    } while (accBytesRead < bufferSize);

    return accBytesRead;
}

int write_to_session_peer(
    SSL*& ssl_session,
    const void* payload,
    size_t payload_length)
{
    int bytes_written = 0;
    int ret = 0;

    while ((bytes_written = SSL_write(ssl_session, payload, payload_length)) <=
        0)
    {
        int error = SSL_get_error(ssl_session, bytes_written);
        if (error == SSL_ERROR_WANT_WRITE)
            continue;
        t_print("Failed! SSL_write returned %d\n", error);
        ret = bytes_written;
        goto exit;
    }

exit:
    return ret;
}
