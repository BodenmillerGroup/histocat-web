"""Add gate table

Revision ID: aba61a88be54
Revises: 56b78600703e
Create Date: 2020-09-30 11:24:30.348903

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = 'aba61a88be54'
down_revision = '56b78600703e'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'gate',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('dataset_id', sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('acquisition_ids', ARRAY(sa.Integer())),
        sa.Column('indices', ARRAY(sa.Integer())),
        sa.Column('cell_ids', ARRAY(sa.Integer())),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('gate')
