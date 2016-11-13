from subprocess import Popen, PIPE


def email_user(attachment: str, email_body: str, email_address: str) -> None:
    """
    sends an email to email address with text contents of email_body and attachment
    attached. Email will come from "ec2-User@<ec2-instance-ip-of-aws-instance>

    :param attachment: the file location of the attachment to append to the email
    :param email_body: text to send in the body of the email
    :param email_address: the address to which the email should be sent"""

    if isinstance(email_body, str):
        email_body = email_body.encode()

    email_args = ['mutt', '-e', 'set content_type="text/html"', '-a', '"%s"' % attachment, '-s',
                  'Remote Process', '--', email_address]
    email_process = Popen(email_args, stdin=PIPE)
    email_process.communicate(email_body)
