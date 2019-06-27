"""Add task table

Revision ID: 639adf01602f
Revises: 485aa67ca586
Create Date: 2019-06-18 16:38:03.448708

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = '639adf01602f'
down_revision = '485aa67ca586'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'task',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('state', sa.String(), nullable=False, index=True),
        sa.Column('name', sa.String(), nullable=False, index=True),
        sa.Column('exitcode', sa.Integer(), index=True),
        sa.Column('time', sa.Interval()),
        sa.Column('memory', sa.Integer()),
        sa.Column('cpu_time', sa.Interval()),
        sa.Column('type', sa.String(), index=True),
        sa.Column('is_collection', sa.Boolean(), index=True),
        sa.Column('parent_id', sa.BigInteger(), index=True),
        sa.Column('data', sa.LargeBinary()),
        sa.Column('submission_id', sa.BigInteger(), sa.ForeignKey("submission.id", ondelete="CASCADE"), index=True),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('task')
