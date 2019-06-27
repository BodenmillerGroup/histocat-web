"""Add slide table

Revision ID: a596fd31b7bd
Revises: 7bed832a38b0
Create Date: 2019-04-10 15:30:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID

# revision identifiers, used by Alembic.

revision = 'a596fd31b7bd'
down_revision = '7bed832a38b0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'slide',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True),
        sa.Column('uid', UUID(), index=True),
        sa.Column('description', sa.String()),
        sa.Column('filename', sa.String(4096)),
        sa.Column('slide_type', sa.String()),
        sa.Column('width_um', sa.Integer()),
        sa.Column('height_um', sa.Integer()),
        sa.Column('image_end_offset', sa.BigInteger()),
        sa.Column('image_start_offset', sa.BigInteger()),
        sa.Column('image_file', sa.String(4096)),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('slide')
